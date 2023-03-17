#include "catsc/cats.h"

#include <memory.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include <float.h>

#ifndef CATS_NPOINTS_REALLOC_STRIDE
#   define CATS_NPOINTS_REALLOC_STRIDE 8
#endif

#ifndef CATCATS_BACKREF_REALLOC_STRIDE
#   define CATS_BACKREF_REALLOC_STRIDE 32
#endif

#ifndef CATSC_WEIGHTS_POOL_STRIDE
#   define CATSC_WEIGHTS_POOL_STRIDE 64
#endif

/*#define COLLECTION_DSTACK_DEBUG 1*/

/* FWD */
struct Cell;

/* A list of references to cells with lazy re-allocation to mitigate heap
 * usage
 * */
struct CellRefs {
    struct Cell ** cells;
    cats_Weight_t * weights;  /* TODO: not needed for cell refs, need for weighted graph for left neighbs */
    size_t nUsed, nAllocated;
};

static void
_cell_refs_init( struct CellRefs * ptr ) {
    assert(ptr);
    ptr->cells = NULL;
    ptr->weights = NULL;
    ptr->nUsed = ptr->nAllocated = 0;
}

/* NOTE: does not (re)allocate weights; in case of reentrant usage, when
 * weights number exceeds previously allocated cells, weights array just gets
 * freed */
static int
_cell_refs_add( struct CellRefs * ptr
              , struct Cell * cellPtr
              ) {
    if( ptr->nUsed + 1 > ptr->nAllocated ) {
        if( SIZE_MAX - CATS_BACKREF_REALLOC_STRIDE <= ptr->nAllocated ) {
            return CATSC_ERROR_ALLOC_FAILURE_BACKREF_COUNTER;
        }
        ptr->nAllocated += CATS_BACKREF_REALLOC_STRIDE;
        struct Cell ** newBackRefs
            = (struct Cell **) realloc( ptr->cells
                                      , ptr->nAllocated * sizeof(struct Cell *) );
        if(!newBackRefs) {
            ptr->nAllocated -= CATS_BACKREF_REALLOC_STRIDE;
            return CATSC_ERROR_ALLOC_FAILURE_BACKREF;
        }
        if(ptr->weights) {
            ptr->weights = (cats_Weight_t *) realloc( ptr->weights,
                    ptr->nAllocated * sizeof(cats_Weight_t *) );
        }
        ptr->cells = newBackRefs;
    }
    ptr->cells[ptr->nUsed] = cellPtr;
    ++(ptr->nUsed);
    return 0;
}

static void
_cell_refs_reset( struct CellRefs * ptr ) {
    ptr->nUsed = 0;
}

static void
_cell_refs_free( struct CellRefs * ptr ) {
    ptr->nUsed = ptr->nAllocated = 0;
    if(ptr->cells) {
        free(ptr->cells);
        ptr->cells = NULL;
    }
    if(ptr->weights) {
        free(ptr->weights);
        ptr->weights = NULL;
    }
}

/*
 * Spatial point
 */

struct cats_Point {
    /** Hit payload data ptr (spatial coordinates) */
    cats_HitData_t data;
    /** Collection of incoming (from left) "cells" */
    struct CellRefs refs;
};

/*
 * CATS Layer representation
 */

struct Layer {
    struct cats_Point * points;
    size_t nPointsAllocated
         , nPointsUsed
         ;
};

static void
_layer_init( struct Layer * l) {
    l->points = NULL;
    l->nPointsUsed = l->nPointsAllocated = 0;
}

static void
_layer_free( struct Layer * l ) {
    if( l->points ) {
        for(size_t nPoint = 0; nPoint < l->nPointsAllocated; ++nPoint) {
            _cell_refs_free( &l->points[nPoint].refs );
        }
        free(l->points);
    }
    l->nPointsUsed = 0;
    l->nPointsAllocated = 0;
}

static int
_layer_add_point( struct Layer * l
                , cats_HitData_t data
                ) {
    if( l->nPointsUsed + 1 > l->nPointsAllocated ) {
        if(SIZE_MAX - CATS_NPOINTS_REALLOC_STRIDE < l->nPointsAllocated) {
            return CATSC_ERROR_ALLOC_FAILURE_POINTS_COUNTER;
        }
        l->nPointsAllocated += CATS_NPOINTS_REALLOC_STRIDE;
        struct cats_Point * newPoints
            = (struct cats_Point *) realloc(l->points, l->nPointsAllocated * sizeof(struct cats_Point));
        if( !newPoints ) {
            /* failed to allocate space */
            l->nPointsAllocated -= CATS_NPOINTS_REALLOC_STRIDE;
            return CATSC_ERROR_ALLOC_FAILURE_POINTS;
        }
        l->points = newPoints;
        /* init point's cell references for newly allocated point */
        for( size_t n = l->nPointsUsed; n < l->nPointsAllocated; ++n ) {
            _cell_refs_init(&l->points[n].refs);
        }
    }

    l->points[l->nPointsUsed].data = data;
    _cell_refs_reset(&l->points[l->nPointsUsed].refs);
    ++(l->nPointsUsed);
    return 0;
}

static void
_layer_reset( struct Layer * l, size_t softLimitHits ) {
    if( softLimitHits && l->nPointsAllocated > softLimitHits ) {
        assert(l->points);
        free(l->points);
        l->points = NULL;
        l->nPointsAllocated = 0;
    }
    l->nPointsUsed = 0;
}


/*
 * CATS layers collection
 */

struct cats_Layers {
    cats_LayerNo_t nLayers;
    struct Layer * layers;
};

/*
 * Public API: layers
 */

struct cats_Layers *
cats_layers_create( cats_LayerNo_t nLayers ) {
    if( nLayers < 3 ) return NULL;  /* nonsense */
    struct cats_Layers * lyrs
        = (struct cats_Layers *) malloc(sizeof(struct cats_Layers));
    if(!lyrs) return NULL;  /* malloc() failed */
    lyrs->layers = (struct Layer *) malloc(nLayers*sizeof(struct Layer));
    for(cats_LayerNo_t n = 0; n < nLayers; ++n) {
        _layer_init(lyrs->layers + n);
    }
    lyrs->nLayers = nLayers;
    return lyrs;
}

void
cats_layers_delete(struct cats_Layers * lyrs) {
    assert(lyrs);
    assert(lyrs->layers);
    for( cats_LayerNo_t n = 0; n < lyrs->nLayers; ++n ) {
        _layer_free(lyrs->layers + n);
    }
    free(lyrs->layers);
    free(lyrs);
}

int
cats_layer_add_point( struct cats_Layers * lsPtr
                    , cats_LayerNo_t layerNo
                    , cats_HitData_t datum
                    ) {
    assert(lsPtr);
    assert(lsPtr->layers);
    if(layerNo >= lsPtr->nLayers) return CATSC_RC_NO_SUCH_LAYER;
    struct Layer * l = lsPtr->layers + layerNo;
    return _layer_add_point( l, datum );
}

size_t
cats_layer_n_points( const struct cats_Layers * lsPtr
                   , cats_LayerNo_t N ) {
    if( !lsPtr ) return 0;
    if( !lsPtr->layers ) return 0;
    if( N >= lsPtr->nLayers ) return 0;
    return lsPtr->layers[N].nPointsUsed;
}

void
cats_layers_reset( struct cats_Layers * ls, size_t softLimitHits ) {
    for( cats_LayerNo_t i = 0; i < ls->nLayers; ++i ) {
        _layer_reset(ls->layers + i, softLimitHits);
    }
}

/*
 * Evaluation
 */

/* A cell is fundamental unit of the algorithm.
 * The structure maintains a pair of hit it refers to and a set of "left"
 * neighbours (from previous layer).
 * */
struct Cell {
    struct cats_Point * from
                    , * to
                    ;
    unsigned int state:31;
    unsigned int doAdvance:1;
    struct CellRefs leftNeighbours;
};

static void
_cell_init( struct Cell * ptr ) {
    ptr->from = NULL;
    ptr->to = NULL;
    ptr->state = 1;
    ptr->doAdvance = 0;
    _cell_refs_init( &ptr->leftNeighbours );
}


struct cats_CellsPool {
    /* Ptr to contiguous block of memory containing all the cells */
    struct Cell * pool;
    size_t nCellsUsed
         , nCellsAllocated;
};

static void
_cell_to_json(struct Cell * cell, FILE * outf) {
    fprintf( outf, "{\"link\":[\"%p\",\"%p\"],\"state\":%u"
           , cell->from, cell->to, cell->state );
    fputs(",\"lneighb\":[", outf);
    for(size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb) {
        if(nNeighb) fputc(',', outf);
        fprintf( outf, "[\"%p\",\"%p\"]"
               , cell->leftNeighbours.cells[nNeighb]->from
               , cell->leftNeighbours.cells[nNeighb]->to
               );
    }
    fputs("]}", outf);
}

static void
_point_to_json( struct cats_Point * p, FILE * outf ) {
    fputs("{\"ptr\":", outf);
    fprintf(outf, "\"%p\",", p);
    fprintf(outf, "\"c\":\"%p\",", p->data);
    fputs("\"refs\":[", outf);
    for( size_t nRef = 0; nRef < p->refs.nUsed; ++nRef ) {
        if(nRef) fputc(',', outf);
        _cell_to_json(p->refs.cells[nRef], outf);
    }
    fputs("]}", outf);
}

void
cats_dump_json( struct cats_Layers * layers
              , FILE * outf ) {
    fputs("[", outf);
    for( cats_LayerNo_t nLayer = 0
       ; nLayer < layers->nLayers
       ; ++nLayer ) {
        if(nLayer) fputs(",", outf);
        fputs("[", outf);
        struct Layer * lPtr = layers->layers + nLayer;
        for( size_t nHit = 0; nHit < lPtr->nPointsUsed; ++nHit ) {
            if(nHit) fputc(',', outf);
            _point_to_json(lPtr->points + nHit, outf);
        }
        fputs("]", outf);
    }
    fputs("]", outf);
}

struct cats_CellsPool *
cats_cells_pool_create(cats_LayerNo_t nLayers) {
    struct cats_CellsPool * a
        = (struct cats_CellsPool *) malloc(sizeof(struct cats_CellsPool));
    if(!a) return NULL;
    a->pool = NULL;
    a->nCellsUsed = a->nCellsAllocated = 0;

    return a;
}

void
cats_cells_pool_reset( struct cats_Layers * ls
                     , struct cats_CellsPool * a
                     , size_t softLimitCells
                     ) {
    for(cats_LayerNo_t nLayer = 0; nLayer < ls->nLayers; ++nLayer) {
        struct Layer * l = ls->layers + nLayer;
        for(size_t nHit = 0; nHit < l->nPointsUsed; ++nHit) {
            struct cats_Point * pt = l->points + nHit;
            _cell_refs_reset(&pt->refs);
        }
    }

    for( size_t nCell = 0; nCell < a->nCellsAllocated; ++nCell ) {
        _cell_refs_reset( &a->pool[nCell].leftNeighbours );
    }

    if( softLimitCells && a->nCellsAllocated > softLimitCells ) {
        /* free anomaluosly large pool set, don't forget about its refs */
        for( size_t nCell = 0; nCell < a->nCellsAllocated; ++nCell ) {
            _cell_refs_free( &a->pool[nCell].leftNeighbours );
        }
        free(a->pool);
        a->pool = NULL;
        a->nCellsAllocated = 0;
    }
    a->nCellsUsed = 0;
}

void
cats_cells_pool_delete( struct cats_CellsPool * a ) {
    if(a->pool) free(a->pool);
    a->nCellsAllocated = a->nCellsUsed = 0;
    free(a);
}

/*                                                          __________________
 * ______________________________________________________ / Graph connections
 */

/* Helper function allocating cells for fully connected graph */
static int
_allocate_fully_connected_cells( struct cats_Layers * ls
                               , struct cats_CellsPool * a
                               , cats_LayerNo_t nMissingLayers
                               , size_t * nCellsNeed_
                               ) {
    /* Estimated number of cells needed */
    size_t nCellsNeed = 0;
    for( cats_LayerNo_t fromNLayer = 0
       ; fromNLayer < ls->nLayers - 1
       ; ++fromNLayer ) {
        for( cats_LayerNo_t toNLayer = fromNLayer + 1
           ; toNLayer < fromNLayer + 2 + nMissingLayers && toNLayer < ls->nLayers
           ; ++toNLayer ) {
            nCellsNeed += ls->layers[fromNLayer].nPointsUsed
                        * ls->layers[toNLayer].nPointsUsed
                        ;
        }
    }
    if(0 == nCellsNeed)
        return CATSC_RC_EMPTY_GRAPH;
    *nCellsNeed_ = nCellsNeed;
    /* (re)allocate/assure required number of cells available */
    if(a->nCellsAllocated < nCellsNeed) {
        cats_cells_pool_reset(ls, a, 0);
        struct Cell * newCells
            = (struct Cell *) realloc( a->pool
                                     , nCellsNeed * sizeof(struct Cell) );
        if( !newCells ) return CATSC_ERROR_ALLOC_FAILURE_CELLS;
        a->pool = newCells;
        /* initialize newly allocated cells for use */
        for( size_t nCell = a->nCellsAllocated
           ; nCell < nCellsNeed
           ; ++nCell ) {
            _cell_init(a->pool + nCell);
        }
        a->nCellsAllocated = nCellsNeed;
    }
    return 0;
}

/* Modifies cells and layer references in a way to reflect the connections
 *
 * - toLayerCellRefs is used as destination collection of to-layer cells refs
 * - "current cell" pointer is used to allocate (borrow) new cells
 * - from/to layers are provided to retrieve collections of hits
 * - test+userdata are for filtering
 * */
static struct Cell *
_connect_layers( struct Cell * cCell
               , struct Layer * fromLayer
               , struct Layer * toLayer
               , cats_Filter_t test_triplet
               , void * userData
               , int * catscErrNo
               ) {
    assert(cCell);
    assert(catscErrNo);
    /* build the cells as cartesian product */
    for(size_t toNHit = 0; toNHit < toLayer->nPointsUsed; ++toNHit) {
        struct cats_Point * toHit = toLayer->points + toNHit;
        for(size_t fromNHit = 0; fromNHit < fromLayer->nPointsUsed; ++fromNHit) {
            struct cats_Point * fromHit = fromLayer->points + fromNHit;
            cCell->from = fromHit;
            cCell->to = toHit;
            cCell->state = 1;
            cCell->doAdvance = 0;
            /* Cells (doublets) in left ("from") layer are created since we
             * iterate from left to right; append the neghbours list. If
             * `fromLayerCellRefs` is null we are probably on the leftmost
             * layer */
            for( size_t i = 0; i < fromHit->refs.nUsed; ++i ) {
                struct Cell * leftCellPtr = fromHit->refs.cells[i];
                assert( leftCellPtr );
                assert( leftCellPtr->to == fromHit );
                /* Check triplet and create cells if test passed */
                if(test_triplet( leftCellPtr->from->data
                               , fromHit->data
                               , toHit->data
                               , userData
                               )) {
                    if(!!(*catscErrNo = _cell_refs_add(&cCell->leftNeighbours, fromHit->refs.cells[i]))) {
                        return NULL;  /* allocation failure */
                    }
                }
            }
            /* Unconditionally append to-hit refs */
            if(!!(*catscErrNo = _cell_refs_add(&toHit->refs, cCell))) {
                return NULL;  /* allocation failure */
            }
            ++cCell;
        }
    }
    return cCell;
}

int
cats_connect( struct cats_Layers * ls
            , struct cats_CellsPool * a
            , cats_Filter_t test_triplet
            , void * userData
            , cats_LayerNo_t nMissingLayers
            ) {
    assert(ls);
    assert(a);
    assert(ls->nLayers > 2);
    ++nMissingLayers;
    size_t nCellsNeed = 0;
    int rc = _allocate_fully_connected_cells(ls, a, nMissingLayers, &nCellsNeed);
    if(rc) return rc;
    /* Create cells */
    struct Cell * cCell = a->pool;  /* tracks last inserted cell */
    for( cats_LayerNo_t toNLayer = 1
       ; toNLayer < ls->nLayers
       ; ++toNLayer ) {
        //cats_LayerNo_t toNLayer = fromNLayer + cNMissing;
        //assert(toNLayer < ls->nLayers);
        for( cats_LayerNo_t fromNLayer = toNLayer - 1
           ; toNLayer - fromNLayer <= nMissingLayers
           ; --fromNLayer ) {
            rc = 0;
            assert(cCell - a->pool < nCellsNeed);  /* pool depleted */
            /* Create cells connecting this layer pair */
            cCell = _connect_layers( cCell
                                   , ls->layers + fromNLayer
                                   , ls->layers + toNLayer
                                   , test_triplet
                                   , userData
                                   , &rc
                                   );
            if(NULL == cCell) {
                assert(rc);  /* _connect_layers() did not return error code */
                return rc;
            }
            if(0 == fromNLayer) break;
        }
    }
    return 0;
}

struct WeightedNeighbour {
    cats_Weight_t weight;
    struct Cell * cell;
};

static int
_compare_weighted_neighbours(const void * a, const void * b) {
    const cats_Weight_t wDiff
         = ( ((const struct WeightedNeighbour *) b)->weight
           - ((const struct WeightedNeighbour *) a)->weight
           );
    /* NOTE: sic! inverting meaning as we want weight to be sorted
     *       descending (first are more reliable ones) */
    if(wDiff > 0) return  1;
    if(wDiff < 0) return -1;
    return 0;
}

/* Modifies cells and layer references in a way to reflect the connections
 *
 * - toLayerCellRefs is used as destination collection of to-layer cells refs
 * - "current cell" pointer is used to allocate (borrow) new cells
 * - from/to layers are provided to retrieve collections of hits
 * - test+userdata are for filtering
 * */
static struct Cell *
_connect_layers_w( struct Cell * cCell
                 , struct Layer * fromLayer
                 , struct Layer * toLayer
                 , cats_WeightedFilter_t test_triplet_w
                 , void * userData
                 , struct WeightedNeighbour ** weightBufferPtr
                 , size_t * weightBufferSizePtr
                 , int * catscErrNo
                 ) {
    assert(cCell);
    assert(catscErrNo);
    /* build the cells as cartesian product
     *  .\\ |     |
     *  .== o === o
     *  .// |  ^  |
     *   /  ^  |  ^- "to hit" on "to layer"
     *   ^  |  `---- "current cell"
     *   |  `------- "from hit" on "from layer"
     *   `---------- "left neghbours" of "current cell"
     * */
    for(size_t toNHit = 0; toNHit < toLayer->nPointsUsed; ++toNHit) {
        struct cats_Point * toHit = toLayer->points + toNHit;
        for(size_t fromNHit = 0; fromNHit < fromLayer->nPointsUsed; ++fromNHit) {
            struct cats_Point * fromHit = fromLayer->points + fromNHit;
            /* init current cell */
            cCell->from = fromHit;
            cCell->to = toHit;
            cCell->state = 1;
            cCell->doAdvance = 0;
            /* number of (weighted) left neighbours for current cell */
            size_t nCurrentWeights = 0;
            /* Cells (doublets) in left ("from") layer are created since we
             * iterate from left to right; append the neghbours list. If
             * `fromLayerCellRefs` is null we are probably on the leftmost
             * layer.
             * Note, that after neighbours are added, they get sorted by
             * weight *within a current tier* (so neighbours sorting are
             * affected by topology -- closes are first). */
            for( size_t i = 0; i < fromHit->refs.nUsed; ++i ) {
                struct Cell * leftCellPtr = fromHit->refs.cells[i];
                assert( leftCellPtr );
                assert( leftCellPtr->to == fromHit );
                /* Check triplet and create cells if test passed */
                cats_Weight_t w = test_triplet_w( leftCellPtr->from->data
                                                , fromHit->data
                                                , toHit->data
                                                , userData
                                                );
                if(w <= 0) continue;  /* no connection for negative weight */
                /* append weights array */
                if(nCurrentWeights == *weightBufferSizePtr) {
                    struct WeightedNeighbour * newBufPtr
                        = (struct WeightedNeighbour *) realloc( *weightBufferPtr
                                                  , sizeof(struct WeightedNeighbour)
                                                        *(*weightBufferSizePtr + CATS_BACKREF_REALLOC_STRIDE)
                                                  );
                    if( !newBufPtr ) {
                        *catscErrNo = CATSC_ERROR_ALLOC_FAILURE_WEIGHTS;
                        return NULL;  /* allocation failure */
                    }
                    *weightBufferPtr = newBufPtr;
                    *weightBufferSizePtr += CATS_BACKREF_REALLOC_STRIDE;
                }
                (*weightBufferPtr)[nCurrentWeights].weight = w;
                (*weightBufferPtr)[nCurrentWeights].cell = fromHit->refs.cells[i];
                ++nCurrentWeights;
            }
            /* Unconditionally append to-hit refs */
            if(!!(*catscErrNo = _cell_refs_add(&toHit->refs, cCell))) {
                return NULL;  /* allocation failure */
            }
            if(0 == nCurrentWeights) {  /* no left neighbours to sort */
                ++cCell;
                continue;
            }
            if(1 == nCurrentWeights) {  /* add single neighbour and that's it */
                if(!!(*catscErrNo = _cell_refs_add(&cCell->leftNeighbours, (*weightBufferPtr)[0].cell))) {
                    return NULL;  /* allocation failure */
                }
                cCell->leftNeighbours.weights
                    = (cats_Weight_t *) malloc(sizeof(cats_Weight_t));  // TODO: pool?
                cCell->leftNeighbours.weights[0] = (*weightBufferPtr)[0].weight;
                ++cCell;
                continue;
            }
            /* sort left neigbours by weights (within the current tier of
             * connectivity) */
            qsort( *weightBufferPtr
                 , nCurrentWeights
                 , sizeof(struct WeightedNeighbour)
                 , _compare_weighted_neighbours );
            cCell->leftNeighbours.weights
                = (cats_Weight_t *) malloc(nCurrentWeights*sizeof(cats_Weight_t));  // TODO: pool?
            assert(cCell->leftNeighbours.weights);
            /* copy left neighbours, in order */
            for(size_t nw = 0; nw < nCurrentWeights; ++nw) {
                if(!!(*catscErrNo = _cell_refs_add(&cCell->leftNeighbours, (*weightBufferPtr)[nw].cell))) {
                    return NULL;  /* allocation failure */
                }
                assert(nw < cCell->leftNeighbours.nUsed);
                assert(cCell->leftNeighbours.weights);
                cCell->leftNeighbours.weights[nw] = (*weightBufferPtr)[nw].weight;
            }
            
            /* shift "current cell" pointer */
            ++cCell;
        }
    }
    return cCell;
}

int
cats_connect_w( struct cats_Layers * ls
              , struct cats_CellsPool * a
              , cats_WeightedFilter_t test_triplet
              , void * userData
              , cats_LayerNo_t nMissingLayers
              ) {
    assert(ls);
    assert(a);
    assert(ls->nLayers > 2);
    ++nMissingLayers;
    size_t nCellsNeed = 0;
    int rc = _allocate_fully_connected_cells(ls, a, nMissingLayers, &nCellsNeed);
    if(rc) return rc;
    /* Create cells */
    struct Cell * cCell = a->pool;  /* tracks last inserted cell */
    size_t weightBufferNItems = 0;
    struct WeightedNeighbour * weightBuffer = NULL;
    for( cats_LayerNo_t toNLayer = 1
       ; toNLayer < ls->nLayers
       ; ++toNLayer ) {
        for( cats_LayerNo_t fromNLayer = toNLayer - 1
           ; toNLayer - fromNLayer <= nMissingLayers
           ; --fromNLayer ) {
            rc = 0;
            assert(cCell - a->pool < nCellsNeed);  /* pool depleted */
            /* Create cells connecting this layer pair */
            cCell = _connect_layers_w( cCell
                                     , ls->layers + fromNLayer
                                     , ls->layers + toNLayer
                                     , test_triplet
                                     , userData
                                     , &weightBuffer
                                     , &weightBufferNItems
                                     , &rc
                                     );
            if(NULL == cCell) {
                assert(rc);  /* _connect_layers() did not return error code */
                if(weightBuffer) free(weightBuffer);
                return rc;
            }
            if(0 == fromNLayer) break;
        }
    }
    if(weightBuffer) free(weightBuffer);
    return 0;
}


static int
_advance_cell_state_if_need(struct Cell * cellPtr) {
    for( size_t nRef = 0
       ; nRef < cellPtr->leftNeighbours.nUsed
       ; ++nRef ) {
        if( cellPtr->leftNeighbours.cells[nRef]->state != cellPtr->state ) continue;
        cellPtr->doAdvance = 1;
    }
    return cellPtr->doAdvance;
}

int
cats_evaluate( struct cats_Layers * ls, FILE * debugJSONStream) {
    assert(ls);
    assert(ls->nLayers > 2);
    /* Evaluate: iterate over cells interfacing the layer 2, 3... and
     * advanced the states of every cell which state matches the state of any
     * of its "left neighbour" */
    int advanced;
    do {
        advanced = 0;
        for( cats_LayerNo_t toNLayer = 0
           ; toNLayer < ls->nLayers
           ; ++toNLayer ) {
            for( size_t nPoint = 0; nPoint < ls->layers[toNLayer].nPointsUsed; ++nPoint ) {
                struct cats_Point * pointPtr = ls->layers[toNLayer].points + nPoint;
                for( size_t nCell = 0; nCell < pointPtr->refs.nUsed; ++nCell ) {
                    struct Cell * cellPtr = pointPtr->refs.cells[nCell];
                    assert( cellPtr->to == pointPtr );
                    advanced |= _advance_cell_state_if_need(cellPtr);
                }
            }
        }
        if(!advanced) {
            if(debugJSONStream) { cats_dump_json(ls, debugJSONStream); }
            break;
        }
        /* update states */
        for( cats_LayerNo_t toNLayer = 1
           ; toNLayer < ls->nLayers
           ; ++toNLayer ) {
            for( size_t nPoint = 0; nPoint < ls->layers[toNLayer].nPointsUsed; ++nPoint ) {
                struct cats_Point * pointPtr = ls->layers[toNLayer].points + nPoint;
                for( size_t nCell = 0; nCell < pointPtr->refs.nUsed; ++nCell ) {
                    struct Cell * cellPtr = pointPtr->refs.cells[nCell];
                    if(cellPtr->doAdvance) {
                        ++(cellPtr->state);
                        cellPtr->doAdvance = 0;
                    }
                }
            }
        }
        if(debugJSONStream) {
            cats_dump_json(ls, debugJSONStream);
            putc(',', debugJSONStream);
        }
    } while(advanced);
    return 0;
}

void
reset_collection_flags( struct cats_CellsPool * pool ) {
    for( size_t n = 0; n < pool->nCellsUsed; ++n ) {
        pool->pool[n].doAdvance = 0x0;
    }
}

/*                                                   _________________________
 * _______________________________________________ / Collection strategies aux
 */

struct PointsStack {
    cats_HitData_t * data;
    ssize_t nTop;
};

#if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
static size_t _gRootLayerNo;
static size_t _gRootNHit;
static size_t _gRootNLink;
static size_t _gRootNMaxLinks;
static size_t _gDbgStack[256][3];

static void _print_stack(size_t N) {
    printf( "ROOT:%zu[%zu] (%zu/%zu)"
          , _gRootLayerNo, _gRootNHit, _gRootNLink, _gRootNMaxLinks );
    for(size_t i = 1; i < N; ++i) {
        printf( " #%zu(%zu/%zu)", _gDbgStack[i][0], _gDbgStack[i][1], _gDbgStack[i][2] );
    }
    fputc('\n', stdout);
}
#endif

void
_stack_push( struct PointsStack * stack
           , cats_HitData_t datum ) {
    stack->data[++(stack->nTop)] = datum;
    assert(stack->nTop >= 0);
}

cats_HitData_t
_stack_top( struct PointsStack * stack ) {
    return stack->data[stack->nTop];
}

cats_HitData_t
_stack_pull( struct PointsStack * stack ) {
    cats_HitData_t datum = _stack_top(stack);
    --(stack->nTop);
    assert(stack->nTop >= -1);
    return datum;
}

/*                                            ________________________________
 * ________________________________________ / Unweighted collection strategies
 */

static void
_eval_from_excessive( struct Cell * cell
                     , struct PointsStack * stack
                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                     , void * userdata
                     , unsigned int minLength
                     ) {
    assert(_stack_top(stack) == cell->to->data);
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state
     && stack->nTop >= minLength ) {
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
    } else {
        for(size_t expectedDiff = 1; expectedDiff <= cell->state; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                if( cell->state - neighb->state == expectedDiff ) {
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    assert( neighb->to == cell->from );
                    _eval_from_excessive(neighb, stack, callback, userdata, minLength);
                }
            }
        }
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

int
cats_visit_dfs_excessive( struct cats_Layers * ls
                        , unsigned int minLength
                        , void (*callback)(const cats_HitData_t *, size_t, void *)
                        , void * userdata
                        ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack;
    stack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = -1;
    /* perform DF starting from the rightmost layer (shall contain highest
     * states, omitting the leftmost one */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_excessive(cell, &stack, callback, userdata, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    return 0;
}

int
cats_visit_dfs_excessive_w( struct cats_Layers * ls
                          , unsigned int minLength
                          , void (*callback)(const cats_HitData_t *, size_t, void *)
                          , void * userdata
                          ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack;
    stack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = -1;

    size_t nNeighbAllocated = CATS_BACKREF_REALLOC_STRIDE;
    struct WeightedNeighbour * wneighb
        = (struct WeightedNeighbour *) malloc(sizeof(struct WeightedNeighbour)*nNeighbAllocated);
    /* perform DF starting from the rightmost layer (shall contain highest
     * states, omitting the leftmost one. For every layer collect cells for
     * sorted left */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        /* For given layer, find highest link state to guarantee sorting
         * within a certain state tier */
        unsigned int maxStateOnThisLayer = 0;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(maxStateOnThisLayer < cell->state) maxStateOnThisLayer = cell->state;
            }
        }
        for(unsigned int cState = maxStateOnThisLayer; cState >= minLength; --cState) {
            assert(cState > 1);
            size_t nNeighb = 0;
            /* prepare array of links of current tier to start with */
            for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
                struct cats_Point * ptStart = l->points + nHit;
                for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                    struct Cell * cell = ptStart->refs.cells[nLink];
                    if(cell->state != cState) continue;
                    assert(cell->state > 1);
                    assert(cell->leftNeighbours.nUsed);  /* state > 1, must have left neighbs */
                    assert(cell->leftNeighbours.weights);  /* graph is not weighted */
                    /* otherwise, copy cell to sort */
                    if(nNeighb == nNeighbAllocated) {
                        struct WeightedNeighbour * newwneighb
                            = (struct WeightedNeighbour *)
                                realloc( wneighb
                                       , sizeof(struct WeightedNeighbour)*(nNeighbAllocated + CATS_BACKREF_REALLOC_STRIDE)
                                       );
                        if(!newwneighb) {
                            return CATSC_ERROR_ALLOC_FAILURE_WEIGHTS;
                        }
                        nNeighbAllocated += CATS_BACKREF_REALLOC_STRIDE;
                        wneighb = newwneighb;
                    }
                    wneighb[nNeighb].cell = cell;
                    wneighb[nNeighb].weight = cell->leftNeighbours.weights[0];
                    ++nNeighb;
                }
            }
            printf(" xxx %zu\n", nNeighb);  // XXX
            if(0 == nNeighb) continue;
            if(nNeighb > 1)
                qsort( wneighb
                     , nNeighb
                     , sizeof(struct WeightedNeighbour)
                     , _compare_weighted_neighbours );
            for(size_t i = 0; i < nNeighb; ++i) {
                struct Cell * cCell = wneighb[i].cell;
                _stack_push(&stack, cCell->to->data);
                _eval_from_excessive(cCell, &stack, callback, userdata, minLength);
                #ifdef NDEBUG
                _stack_pull(&stack);
                #else
                assert(_stack_pull(&stack) == cCell->to->data);
                #endif
            }
        }
    }
    free(stack.data);
    free(wneighb);
    return 0;
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from( struct Cell * cell
          , struct PointsStack * stack
          , void (*callback)(const cats_HitData_t *, size_t, void *)
          , void * userdata
          , unsigned int minLength
          ) {
    assert(_stack_top(stack) == cell->to->data);
    cell->doAdvance = 1;
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
    } else {
        for(size_t expectedDiff = 1; expectedDiff < cell->state; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                //if( neighb->doAdvance ) continue;  // TODO?
                if( cell->state - neighb->state == expectedDiff ) {
                    assert( neighb->to == cell->from );
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    _eval_from(neighb, stack, callback, userdata, minLength);
                }
            }
        }
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

int
cats_visit_dfs_moderate( struct cats_Layers * ls
                       , unsigned int minLength
                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                       , void * userdata
                       ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack;
    stack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = -1;
    /* evaluate traversing */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from(cell, &stack, callback, userdata, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    return 0;
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from_strict( struct Cell * cell
                 , struct PointsStack * stack
                 , struct PointsStack * prevStack
                 , void (*callback)(const cats_HitData_t *, size_t, void *)
                 , void * userdata
                 , unsigned int minLength
                 ) {
    assert(_stack_top(stack) == cell->to->data);
    cell->doAdvance = 1;
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        /* Compare current stack with previous successful one:
         *  - if current is longer or has equal length, commit it (must be
         *    different because of the graph properties)
         *  - if current has at least one distinct element from prevous, commit
         *    it. */
        if( prevStack->nTop <= stack->nTop || stack->nTop < 0 ) goto commit;
        for( ssize_t i = 0; i <= stack->nTop; ++i ) {
            int found = 0;
            for( ssize_t j = 0; j <= prevStack->nTop; ++j ) {
                if( prevStack->data[j] != stack->data[i] ) continue;
                found = 1;
                break;
            }
            if(found) continue;
            goto commit;  /* current stack brought at least one uniq element */
        }
        goto omit_subsequence;
commit:
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        /* Copy stack for comparison */
        prevStack->nTop = stack->nTop;
        for(ssize_t n = 0; n <= stack->nTop; ++n) prevStack->data[n] = stack->data[n];
    } else {
        for(size_t expectedDiff = 1; expectedDiff < cell->state; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                if( cell->state - neighb->state == expectedDiff ) {
                    assert( neighb->to == cell->from );
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    _eval_from_strict(neighb, stack, prevStack, callback, userdata, minLength);
                }
            }
        }
    }
omit_subsequence:
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

int
cats_visit_dfs_strict( struct cats_Layers * ls
                     , unsigned int minLength
                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                     , void * userdata
                     ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack
                     , prevStack
                     ;
    stack.data =     (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    prevStack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = prevStack.nTop = -1;
    /* evaluate traversing */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_strict(cell, &stack, &prevStack, callback, userdata, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    free(prevStack.data);
    return 0;
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from_l( struct Cell * cell
            , struct PointsStack * stack
            , struct PointsStack * prevStack
            , void (*callback)(const cats_HitData_t *, size_t, void *)
            , void * userdata
            , unsigned int nMissingLayers
            , unsigned int minLength
            ) {
    assert(_stack_top(stack) == cell->to->data);
    cell->doAdvance = 1;  /* mark as visited ("remove" from set) */
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        /* same as for the strict */
        if( prevStack->nTop <= stack->nTop || stack->nTop < 0 ) goto commit;
        for( ssize_t i = 0; i <= stack->nTop; ++i ) {
            int found = 0;
            for( ssize_t j = 0; j <= prevStack->nTop; ++j ) {
                if( prevStack->data[j] != stack->data[i] ) continue;
                found = 1;
                break;
            }
            if(found) continue;
            goto commit;  /* current stack brought at least one uniq element */
        }
        goto omit_subsequence;
commit:
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        /* Copy stack for comparison */
        prevStack->nTop = stack->nTop;
        for(ssize_t n = 0; n <= stack->nTop; ++n) prevStack->data[n] = stack->data[n];
    } else {
        for(size_t expectedDiff = 1; expectedDiff < cell->state /*nMissingLayers + 1*/; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                if( neighb->doAdvance ) continue;  /* omit visited ("removed") */
                if( cell->state - neighb->state == expectedDiff ) {
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    assert( neighb->to == cell->from );
                    _eval_from_l(neighb, stack, prevStack, callback, userdata, nMissingLayers, minLength);
                }
            }
        }
    }
    //if( stack->nTop >= minLength ) {
    //    callback(stack->data, stack->nTop + 1, userdata);
    //    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
    //    _print_stack(stack->nTop + 1);
    //    #endif
    //}
omit_subsequence:
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

int
cats_visit_dfs_longest( struct cats_Layers * ls
                      , unsigned int minLength
                      , unsigned int nMissingLayers
                      , void (*callback)(const cats_HitData_t *, size_t, void *)
                      , void * userdata
                      ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack
                     , prevStack
                     ;
    stack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    prevStack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = prevStack.nTop = -1;
    /* evaluate traversing */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // omit "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_l(cell, &stack, &prevStack, callback, userdata, nMissingLayers, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    free(prevStack.data);
    return 0;
}

/*                          * * *   * * *   * * *                            */

static int
_eval_from_winning_l( struct Cell * cell
                    , struct PointsStack * stack
                    , void (*callback)(const cats_HitData_t *, size_t, void *)
                    , void * userdata
                    , unsigned int nMissingLayers
                    , unsigned int minLength
                    ) {
    assert(_stack_top(stack) == cell->to->data);
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        #ifdef NDEBUG
        _stack_pull(stack);
        #else
        assert(_stack_pull(stack) == cell->from->data);
        #endif
        cell->doAdvance = 1;  /* mark as visited ("remove" from set) */
        return 1;
    } else {
        for(size_t expectedDiff = 1; expectedDiff <= nMissingLayers + 1; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                if( neighb->doAdvance ) continue;  /* omit visited ("removed") */
                if( cell->state - neighb->state == expectedDiff ) {
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    assert( neighb->to == cell->from );
                    int won = _eval_from_winning_l(neighb, stack, callback,
                                userdata, nMissingLayers, minLength);
                    if(won) {
                        cell->doAdvance = 1;  /* mark as visited ("remove" from set) */
                        #ifdef NDEBUG
                        _stack_pull(stack);
                        #else
                        assert(_stack_pull(stack) == cell->from->data);
                        #endif
                        return won;
                    }
                }
            }
        }
    }
    if( stack->nTop >= minLength ) {
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        #ifdef NDEBUG
        _stack_pull(stack);
        #else
        assert(_stack_pull(stack) == cell->from->data);
        #endif
        cell->doAdvance = 1;  /* mark as visited ("remove" from set) */
        return 1;
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
    return 0;
}

int
cats_visit_dfs_winning( struct cats_Layers * ls
                      , unsigned int minLength
                      , unsigned int nMissingLayers
                      , void (*callback)(const cats_HitData_t *, size_t, void *)
                      , void * userdata
                      ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack;
    stack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = -1;
    /* evaluate traversing */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                if(cell->doAdvance) continue;  // omit "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                /*int won =*/ _eval_from_winning_l( cell, &stack, callback
                        , userdata, nMissingLayers, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    return 0;
}

/* ... TODO: unweighted BFS visitors */

/*                                              ______________________________
 * __________________________________________ / Weighted collection strategies
 */

/* ... TODO: weighted wcessiva and moderate implems */

/*                          * * *   * * *   * * *                            */

#if 0
static void
_eval_from_strict_w( struct Cell * cell
                   , struct PointsStack * stack
                   , struct PointsStack * prevStack
                   , void (*callback)(const cats_HitData_t *, size_t, void *)
                   , void * userdata
                   , cats_WeightedFilter_t wFilter
                   , void * wFilterUserdata
                   , unsigned int minLength
                 ) {
    assert(_stack_top(stack) == cell->to->data);
    cell->doAdvance = 1;
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        /* Compare current stack with previous successful one:
         *  - if current is longer or has equal length, commit it (must be
         *    different because of the graph properties)
         *  - if current has at least one distinct element from prevous, commit
         *    it. */
        if( prevStack->nTop <= stack->nTop || stack->nTop < 0 ) goto commit;
        for( ssize_t i = 0; i <= stack->nTop; ++i ) {
            int found = 0;
            for( ssize_t j = 0; j <= prevStack->nTop; ++j ) {
                if( prevStack->data[j] != stack->data[i] ) continue;
                found = 1;
                break;
            }
            if(found) continue;
            goto commit;  /* current stack brought at least one uniq element */
        }
        goto omit_subsequence;
commit:
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        /* Copy stack for comparison */
        prevStack->nTop = stack->nTop;
        for(ssize_t n = 0; n <= stack->nTop; ++n) prevStack->data[n] = stack->data[n];
    } else {
        if( ! cell->leftNeighbours.weights ) {
            /* Important note: to implement strict strategy based on previous
             * stack we should preserve a feature of connected cells: closer
             * cells are coming first in the "leftNeighbours" list. Otherwise,
             * insertion sequence get broken and, say (3,1) (2,1) (1,1) with weight
             * 0.1 will be visited (3,1) (1,1) with weight 0.9. To perform this we
             * sort connections by tiers of their connectivity. */
            cell->leftNeighbours.weights 
                = (cats_Weight_t *) malloc(sizeof(cats_Weight_t)*cell->leftNeighbours.nUsed);
            // ...
        }   
        for(size_t expectedDiff = 1; expectedDiff < cell->state; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                if( cell->state - neighb->state == expectedDiff ) {
                    assert( neighb->to == cell->from );
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                    _eval_from_strict_w( neighb, stack, prevStack
                                       , callback, userdata
                                       , wFilter, wFilterUserdata
                                       , minLength );
                }
            }
        }
    }
omit_subsequence:
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

int
cats_visit_dfs_strict_w( struct cats_Layers * ls
                       , unsigned int minLength
                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                       , void * userdata
                       , cats_WeightedFilter_t wFilter
                       , void * wFilterUserdata
                       ) {
    if(minLength < 3) return CATSC_RC_BAD_MIN_LENGTH;
    --minLength;
    /* initialize traversing stack */
    struct PointsStack stack
                     , prevStack
                     ;
    stack.data =     (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    prevStack.data = (const void **) malloc(ls->nLayers*sizeof(cats_HitData_t));
    stack.nTop = prevStack.nTop = -1;
    /* evaluate traversing */
    for( cats_LayerNo_t nLayer = ls->nLayers - 1; nLayer > 0; --nLayer ) {
        struct Layer * l = ls->layers + nLayer;
        for( size_t nHit = 0; nHit < l->nPointsUsed; ++nHit ) {
            struct cats_Point * ptStart = l->points + nHit;
            _stack_push(&stack, ptStart->data);
            for( size_t nLink = 0; nLink < ptStart->refs.nUsed; ++nLink ) {
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_strict_w( cell, &stack, &prevStack
                                   , callback, userdata
                                   , wFilter, wFilterUserdata
                                   , minLength );
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
    free(prevStack.data);
    return 0;
}
#endif

/*                          * * *   * * *   * * *                            */

