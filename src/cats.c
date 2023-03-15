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
    size_t nUsed, nAllocated;
};

static void
_cell_refs_init( struct CellRefs * ptr ) {
    assert(ptr);
    ptr->cells = NULL;
    ptr->nUsed = ptr->nAllocated = 0;
}

static int
_cell_refs_add( struct CellRefs * ptr, struct Cell * cellPtr ) {
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
    for(size_t fromNHit = 0; fromNHit < fromLayer->nPointsUsed; ++fromNHit) {
        struct cats_Point * fromHit = fromLayer->points + fromNHit;
        for(size_t toNHit = 0; toNHit < toLayer->nPointsUsed; ++toNHit) {
            struct cats_Point * toHit = toLayer->points + toNHit;
            cCell->from = fromHit;
            cCell->to = toHit;
            cCell->state = 1;
            cCell->doAdvance = 0;
            /* Cells (doublets) in left ("from") layer are created since we
             * iterate from left to right; fill up the neghbours list. If
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
cats_evaluate( struct cats_Layers * ls
             , struct cats_CellsPool * a
             , cats_Filter_t test_triplet
             , void * userData
             , unsigned int nMissingLayers
             , FILE * debugJSONStream
             ) {
    assert(ls);
    assert(a);
    assert(ls->nLayers > 2);
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
    /* Create cells */
    struct Cell * cCell = a->pool;  /* tracks last inserted cell */
    #if 1
    #if 1  //assert(nMissingLayers < ls->nLayers);
    for( cats_LayerNo_t toNLayer = 1
       ; toNLayer < ls->nLayers
       ; ++toNLayer ) {
        //cats_LayerNo_t toNLayer = fromNLayer + cNMissing;
        //assert(toNLayer < ls->nLayers);
        for( cats_LayerNo_t fromNLayer = toNLayer - 1
           ; toNLayer - fromNLayer <= nMissingLayers
           ; --fromNLayer ) {
            int rc = 0;
            assert(cCell - a->pool < nCellsNeed);  /* pool depleted */
            printf(" xxx connecting %d -> %d\n", (int) fromNLayer, (int) toNLayer);
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
    #else
    for( cats_LayerNo_t cNMissing = 1; cNMissing <= nMissingLayers; ++ cNMissing ) {
        for( cats_LayerNo_t fromNLayer = 0
           ; fromNLayer < ls->nLayers - cNMissing
           ; ++fromNLayer ) {
            //cats_LayerNo_t toNLayer = fromNLayer + cNMissing;
            //assert(toNLayer < ls->nLayers);
            for( cats_LayerNo_t toNLayer = fromNLayer + 1
               ; toNLayer <= fromNLayer + cNMissing //toNLayer < ls->nLayers
               ; ++toNLayer ) {
                int rc = 0;
                assert(cCell - a->pool < nCellsNeed);  /* pool depleted */
                printf(" xxx connecting %d -> %d\n", (int) fromNLayer, (int) toNLayer);
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
            }
        }
    }
    #endif
    #else
    for( cats_LayerNo_t fromNLayer = 0
       ; fromNLayer < ls->nLayers - 1
       ; ++fromNLayer ) {
        for( cats_LayerNo_t toNLayer = fromNLayer + 1
           ; toNLayer < fromNLayer + 2 + nMissingLayers && toNLayer < ls->nLayers
           ; ++toNLayer ) {
            int rc = 0;
            /* Create cells connecting this layer pair */
            cCell = _connect_layers( cCell
                                   , ls->layers + fromNLayer
                                   , ls->layers + toNLayer
                                   , test_triplet
                                   , userData
                                   , &rc
                                   );
            if(NULL == cCell) {
                assert(rc);
                return rc;
            }
        }
    }
    #endif
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

/*                          * * *   * * *   * * *                            */

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

/*                          * * *   * * *   * * *                            */

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
                printf( "xxx exp.diff=%d <= %d, df=%d\n"
                      , (int) expectedDiff
                      , (int) cell->state
                      , (int) cell->state - neighb->state
                      );  // XXX
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
            puts("next neighb\n");  // XXX
        }
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

void
cats_for_each_track_candidate_excessive( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers  // xxx
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdata
                                       ) {
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
}

void
reset_collection_flags( struct cats_CellsPool * pool ) {
    for( size_t n = 0; n < pool->nCellsUsed; ++n ) {
        pool->pool[n].doAdvance = 0x0;
    }
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

void
cats_for_each_track_candidate( struct cats_Layers * ls
                             , unsigned int minLength
                             , void (*callback)(const cats_HitData_t *, size_t, void *)
                             , void * userdata
                             ) {
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
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from_l( struct Cell * cell
            , struct PointsStack * stack
            , void (*callback)(const cats_HitData_t *, size_t, void *)
            , void * userdata
            , unsigned int nMissingLayers
            , unsigned int minLength
            ) {
    assert(_stack_top(stack) == cell->to->data);
    cell->doAdvance = 1;  /* mark as visited ("remove" from set) */
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state && stack->nTop >= minLength ) {
        callback(stack->data, stack->nTop + 1, userdata);
        #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
        _print_stack(stack->nTop + 1);
        #endif
        /* testing: mark when pushing, mark incoming cells with state <
         * current as visited */
        //for( size_t n = 0; n < cell->to->refs.nUsed; ++n ) {
        //    if( cell->to->refs.cells[n]->state < cell->state )
        //        cell->to->refs.cells[n]->doAdvance = 1;
        //}
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
                    _eval_from_l(neighb, stack, callback, userdata, nMissingLayers, minLength);
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
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

void
cats_for_each_longest_track_candidate( struct cats_Layers * ls
                                     , unsigned int minLength
                                     , unsigned int nMissingLayers
                                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                                     , void * userdata
                                     ) {
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
                if(cell->doAdvance) continue;  // omit "visited"
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_l(cell, &stack, callback, userdata, nMissingLayers, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
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

void
cats_for_each_winning_track_candidate( struct cats_Layers * ls
                                     , unsigned int minLength
                                     , unsigned int nMissingLayers
                                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                                     , void * userdata
                                     ) {
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
                printf( "xxx %d,%d --%p--> : visited=%d\n"                          // XXX
                      , (int) nLayer, (int) nHit, cell, (int) cell->doAdvance);   // XXX
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
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from_w( struct Cell * cell
            , struct PointsStack * stack
            , cats_WeightedFilter_t wf
            , void * userdataFilter
            , void (*callback)(const cats_HitData_t *, size_t, void *)
            , void * userdata
            , unsigned int nMissingLayers
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
        /* Find out best link among concurrent ones for all depth */
        double bestWeight = DBL_MIN;
        struct Cell * bestCell = NULL;
        for(size_t expectedDiff = 1; expectedDiff <= nMissingLayers + 1; ++expectedDiff ) {
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                assert(stack->nTop > -1);
                double w = wf( stack->data[stack->nTop]
                             , neighb->from->data
                             , neighb->to->data
                             , userdataFilter );
                if( w > bestWeight ) {
                    bestWeight = w;
                    bestCell = neighb;
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                }
            }
        }
        /* Proceed with winner */
        if(bestCell) {
            assert( bestCell->to == cell->from );
            _eval_from_w(bestCell, stack, wf, userdataFilter, callback, userdata
                    , nMissingLayers, minLength);
        }
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

void
cats_for_each_track_candidate_w( struct cats_Layers * ls
                               , unsigned int minLength
                               , unsigned int nMissingLayers
                               , cats_WeightedFilter_t wf
                               , void * wfData
                               , void (*callback)(const cats_HitData_t *, size_t, void *)
                               , void * userdataCollect
                               ) {
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
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_w(cell, &stack, wf, wfData, callback, userdataCollect, nMissingLayers, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
}

/*                          * * *   * * *   * * *                            */

static void
_eval_from_w_longest( struct Cell * cell
                    , struct PointsStack * stack
                    , cats_WeightedFilter_t wf
                    , void * userdataFilter
                    , void (*callback)(const cats_HitData_t *, size_t, void *)
                    , void * userdata
                    , unsigned int nMissingLayers
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
        /* Find out best link among concurrent ones for all depth */
        for(size_t expectedDiff = 1; expectedDiff <= nMissingLayers + 1; ++expectedDiff ) {
            double bestWeight = DBL_MIN;
            struct Cell * bestCell = NULL;
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                assert(stack->nTop > -1);
                double w = wf( stack->data[stack->nTop]
                             , neighb->from->data
                             , neighb->to->data
                             , userdataFilter );
                if( w > bestWeight ) {
                    bestWeight = w;
                    bestCell = neighb;
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                }
            }
            /* Proceed with winner at current depth */
            if(bestCell) {
                assert( bestCell->to == cell->from );
                _eval_from_w_longest(bestCell, stack, wf, userdataFilter
                        , callback, userdata, nMissingLayers, minLength);
            }
        }
    }
    #ifdef NDEBUG
    _stack_pull(stack);
    #else
    assert(_stack_pull(stack) == cell->from->data);
    #endif
}

void
cats_for_each_longest_track_candidate_w( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers
                                       , cats_WeightedFilter_t wf
                                       , void * wfData
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdataCollect
                                       ) {
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
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from_w_longest(cell, &stack, wf, wfData, callback
                        , userdataCollect, nMissingLayers, minLength);
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
}

/*                          * * *   * * *   * * *                            */

static int
_eval_from_w_winning( struct Cell * cell
                    , struct PointsStack * stack
                    , cats_WeightedFilter_t wf
                    , void * userdataFilter
                    , void (*callback)(const cats_HitData_t *, size_t, void *)
                    , void * userdata
                    , unsigned int nMissingLayers
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
        #ifdef NDEBUG
        _stack_pull(stack);
        #else
        assert(_stack_pull(stack) == cell->from->data);
        #endif
        return 1;
    } else {
        /* Find out best link among concurrent ones for all depth */
        for(size_t expectedDiff = 1; expectedDiff <= nMissingLayers + 1; ++expectedDiff ) {
            double bestWeight = DBL_MIN;
            struct Cell * bestCell = NULL;
            for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
                struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
                assert(stack->nTop > -1);
                double w = wf( stack->data[stack->nTop]
                             , neighb->from->data
                             , neighb->to->data
                             , userdataFilter );
                if( w > bestWeight ) {
                    bestWeight = w;
                    bestCell = neighb;
                    #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                    _gDbgStack[stack->nTop][0] = expectedDiff;
                    _gDbgStack[stack->nTop][1] = nNeighb;
                    _gDbgStack[stack->nTop][2] = cell->leftNeighbours.nUsed;
                    #endif
                }
            }
            /* Proceed with winner at current depth */
            if(bestCell) {
                assert( bestCell->to == cell->from );
                int won = _eval_from_w_winning( bestCell, stack, wf, userdataFilter
                                              , callback, userdata, nMissingLayers, minLength);
                if(won) {
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

void
cats_for_each_winning_track_candidate_w( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers
                                       , cats_WeightedFilter_t wf
                                       , void * wfData
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdataCollect
                                       ) {
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
                #if defined(COLLECTION_DSTACK_DEBUG) && COLLECTION_DSTACK_DEBUG
                _gRootLayerNo = nLayer;
                _gRootNHit = nHit;
                _gRootNLink = nLink;
                _gRootNMaxLinks = ptStart->refs.nUsed;
                #endif
                struct Cell * cell = ptStart->refs.cells[nLink];
                if(cell->state < minLength) continue;
                // use `doAdvance' flag to mark visited cells
                if(cell->doAdvance) continue;  // "visited"
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                int won = _eval_from_w_winning( cell, &stack, wf, wfData, callback
                                              , userdataCollect, nMissingLayers, minLength);
                if(won) break;
            }
            #ifdef NDEBUG
            _stack_pull(&stack);
            #else
            assert(_stack_pull(&stack) == ptStart->data);
            #endif
        }
    }
    free(stack.data);
}

