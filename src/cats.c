#include "catsc/cats.h"

#include <memory.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#ifndef CATS_NPOINTS_REALLOC_STRIDE
#   define CATS_NPOINTS_REALLOC_STRIDE 8
#endif

#ifndef CATCATS_BACKREF_REALLOC_STRIDE
#   define CATS_BACKREF_REALLOC_STRIDE 8
#endif

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
        ptr->nAllocated += CATS_BACKREF_REALLOC_STRIDE;
        struct Cell ** newBackRefs
            = (struct Cell **) reallocarray( ptr->cells
                                           , ptr->nAllocated
                                           , sizeof(struct Cell *) );
        if(!newBackRefs) {
            ptr->nAllocated -= CATS_BACKREF_REALLOC_STRIDE;
            return -1;
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
    /** Collection of incoming "cells" */
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
        l->nPointsAllocated += CATS_NPOINTS_REALLOC_STRIDE;
        struct cats_Point * newPoints
            = (struct cats_Point *) reallocarray(l->points, l->nPointsAllocated, sizeof(struct cats_Point));
        if( !newPoints ) {
            /* failed to allocate space */
            l->nPointsAllocated -= CATS_NPOINTS_REALLOC_STRIDE;
            return -2;
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
    if(layerNo >= lsPtr->nLayers) return -1;
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
    unsigned int state;
    struct CellRefs leftNeighbours;
};

static void
_cell_init( struct Cell * ptr ) {
    ptr->from = NULL;
    ptr->to = NULL;
    ptr->state = 1;
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
    const size_t nDoubletsGroups = nLayers - 1;
    struct cats_CellsPool * a
        = (struct cats_CellsPool *) malloc(sizeof(struct cats_CellsPool));

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
               ) {
    /* build the cells as cartesian product */
    for(size_t fromNHit = 0; fromNHit < fromLayer->nPointsUsed; ++fromNHit) {
        struct cats_Point * fromHit = fromLayer->points + fromNHit;
        for(size_t toNHit = 0; toNHit < toLayer->nPointsUsed; ++toNHit) {
            struct cats_Point * toHit = toLayer->points + toNHit;
            cCell->from = fromHit;
            cCell->to = toHit;
            cCell->state = 1;
            /* Cells (doublets) in left ("from") layer are created since we
             * iterate from left to right; fill up the neghbours list. If
             * `fromLayerCellRefs` we are probably on the leftmost layer (so
             * current cells have no "left neghbours)" */
            for( size_t i = 0; i < fromHit->refs.nUsed; ++i ) {
                struct Cell * leftCellPtr = fromHit->refs.cells[i];
                assert( leftCellPtr->to == fromHit );
                /* Check triplet and create cells if test passed */
                if(test_triplet( leftCellPtr->from->data
                               , fromHit->data
                               , toHit->data
                               , userData
                               )) {
                    _cell_refs_add(&cCell->leftNeighbours, fromHit->refs.cells[i]);
                }
            }
            /* Unconditionally append to-hit refs */
            _cell_refs_add(&toHit->refs, cCell);
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
        if( cellPtr->leftNeighbours.cells[nRef]->state == cellPtr->state ) {
            ++(cellPtr->state);
            return 1;
        }
    }
    return 0;
}

int
cats_evolve( struct cats_Layers * ls
           , struct cats_CellsPool * a
           , cats_Filter_t test_triplet
           , void * userData
           , unsigned int nMissingLayers ) {
    assert(ls);
    assert(a);
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
    /* (re)allocate/assure required number of cells available */
    if(a->nCellsAllocated < nCellsNeed) {
        cats_cells_pool_reset(ls, a, 0);
        struct Cell * newCells
            = (struct Cell *) reallocarray( a->pool
                                          , nCellsNeed
                                          , sizeof(struct Cell) );
        if( !newCells ) return -1;
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
    struct Cell * cCell = a->pool;
    for( cats_LayerNo_t fromNLayer = 0
       ; fromNLayer < ls->nLayers - 1
       ; ++fromNLayer ) {
        for( cats_LayerNo_t toNLayer = fromNLayer + 1
           ; toNLayer < fromNLayer + 2 + nMissingLayers && toNLayer < ls->nLayers
           ; ++toNLayer ) {
            /* Create cells connecting this layer pair */
            cCell = _connect_layers( cCell
                                   , ls->layers + fromNLayer
                                   , ls->layers + toNLayer
                                   , test_triplet
                                   , userData
                                   );
        }
    }
    /* Evaluate: iterate over cells interfacing the layer 2, 3... and
     * advanced the states of every cell which state matches the state of any
     * of its "left neighbour" */
    int advanced;
    do {
        advanced = 0;
        for( cats_LayerNo_t toNLayer = 1
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
    } while(advanced);
    return 0;
}

struct PointsStack {
    cats_HitData_t * data;
    ssize_t nTop;
};

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

static void
_eval_from( struct Cell * cell
          , struct PointsStack * stack
          , void (*callback)(const cats_HitData_t *, size_t, void *)
          , void * userdata ) {
    assert(_stack_top(stack) == cell->to->data);
    _stack_push(stack, cell->from->data);
    if( 1 == cell->state ) {
        callback(stack->data, stack->nTop + 1, userdata);
    } else {
        for( size_t nNeighb = 0; nNeighb < cell->leftNeighbours.nUsed; ++nNeighb ) {
            struct Cell * neighb = cell->leftNeighbours.cells[nNeighb];
            if( neighb->state + 1 == cell->state ) {
                assert( neighb->to == cell->from );
                _eval_from(neighb, stack, callback, userdata);
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
    /* initialize traversing stack */
    struct PointsStack stack;
    stack.data = malloc(ls->nLayers*sizeof(cats_HitData_t));
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
                assert(cell->leftNeighbours.nUsed); /* cell can not have state >1 without neghbours */
                _eval_from(cell, &stack, callback, userdata);
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

