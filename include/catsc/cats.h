#ifndef H_CATS_H
#define H_CATS_H

#include <stdio.h>

#include "catsc/config.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CATS_HIT_ID_TYPE
#   define CATS_HIT_ID_TYPE unsigned int
#endif

/** Type identifying a hit */
typedef CATS_HIT_ID_TYPE cats_HitID_t;

/** Filter callback type */
typedef int (*cats_Filter_t)( cats_HitID_t, const void *
                            , cats_HitID_t, const void *
                            , cats_HitID_t, const void *
                            , void * );

/** Type to identify layers (internal) */
typedef unsigned int cats_LayerNo_t;

/* FWD */
struct cats_Layers;
struct cats_CellsPool;

/**\brief Create layers collection
 *
 * Requires knowledge of foreseen number of layers in advance */
struct cats_Layers * cats_layers_create( cats_LayerNo_t );

/**\brief Delete automaton */
void cats_layers_delete(struct cats_Layers *);

/**\brief Adds a hit to certain layer
 *
 * Coordinates triplet may not be in agreement with "sorting parameter" choosen
 * for layer.
 *
 * \returns non-zero if failed to allocate memory for new point.
 * */
int cats_layer_add_point( struct cats_Layers *, cats_LayerNo_t
                        , const void *
                        , cats_HitID_t id
                        );

/**\brief Erases the hits collected within layers
 *
 * This will effectively "delete" hits currently being stored in layers
 * collection, making created topology available for reentrant usage.
 * */
void cats_layers_reset(struct cats_Layers *, size_t softLimitHits);

/**\brief Constructor of cells pool
 *
 * Sets up reentrant structure to be used with certain amount of layers.
 * To destroy object use `cats_cells_pool_delete()` */
struct cats_CellsPool * cats_cells_pool_create(cats_LayerNo_t);

/**\brief Dumps state of the automaton as columnar file
 * */
void cats_dump_json( struct cats_Layers *
                   , FILE * );

/**\brief Deletes the automaton.
 *
 * Frees memory occupied during automaton lifecycle.
 * */
void cats_cells_pool_delete( struct cats_CellsPool * );

/**\brief Re-sets the pool for rentrant usage*/
void cats_cells_pool_reset( struct cats_Layers *
                          , struct cats_CellsPool *
                          , size_t softLimitCells
                          );

/**\brief Performs forward evaluation (1st stage) of CATS algorithm
 *
 * Constructs cellular automaton based on hits info and evaluates states till
 * convergence.
 * */
int cats_evolve( struct cats_Layers *
               , struct cats_CellsPool *
               , cats_Filter_t test_triplet
               , void * userData
               , unsigned int nMissingLayers
               );

void
cats_for_each_track_candidate( struct cats_Layers * ls
                             , unsigned int minLength
                             , void (*callback)(cats_HitID_t *, size_t, void *)
                             , void * userdata
                             );

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  /* H_CATS_H */
