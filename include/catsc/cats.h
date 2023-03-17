#ifndef H_CATS_H
#define H_CATS_H

#include <stdio.h>

#include "catsc/config.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CATSC_ERROR_ALLOC_FAILURE_BACKREF_COUNTER    0x1  /* backref counter exceeded */
#define CATSC_ERROR_ALLOC_FAILURE_BACKREF            0x2  /* backref instance allocation failure */
#define CATSC_ERROR_ALLOC_FAILURE_POINTS_COUNTER     0x3  /* pointer counter on layer exceeded */
#define CATSC_ERROR_ALLOC_FAILURE_POINTS             0x4  /* pointer instance allocation failed */
#define CATSC_ERROR_ALLOC_FAILURE_CELLS              0x5  /* failed to (re)allocate cells */
#define CATSC_ERROR_ALLOC_FAILURE_WEIGHTS            0x6  /* failed to (re)allocate weights buffer on layer */
#define CATSC_ERROR_RUNTIME_LOGIC                   0x10  /* runtime error flag */
#define CATSC_RC_NO_SUCH_LAYER                      0x11  /* No layer of given number */
#define CATSC_RC_EMPTY_GRAPH                        0x12  /* Graph has no connections (not actually erorr) */
#define CATSC_RC_BAD_LAYERS_NO                      0x13  /* Layers number is too small */
#define CATSC_RC_BAD_MIN_LENGTH                     0x14  /* Too short min length for collect() */

#ifndef CATS_HIT_DATA_TYPE
#   define CATS_HIT_DATA_TYPE const void*
#endif

#ifndef CATS_WEIGHT_DATA_TYPE
#   define CATS_WEIGHT_DATA_TYPE float
#endif

/** Type of hit data payload */
typedef CATS_HIT_DATA_TYPE cats_HitData_t;

/**\brief Filter callback type
 *
 * This callback shall return non-zero code if given hits triplet can be a
 * valid track piece */
typedef int (*cats_Filter_t)( cats_HitData_t
                            , cats_HitData_t
                            , cats_HitData_t
                            , void * );

typedef CATS_WEIGHT_DATA_TYPE cats_Weight_t;

/**\brief Weighted filter callback type
 *
 * Shall return relative weight for the triplet to be combined in a track
 * piece. Returns zero or negavie value if combination should not be
 * considered.
 * */
typedef cats_Weight_t (*cats_WeightedFilter_t)( cats_HitData_t
                                              , cats_HitData_t
                                              , cats_HitData_t
                                              , void * );

/** Type to identify layers (internal) */
typedef unsigned int cats_LayerNo_t;

/* FWD */
struct cats_Layers;
struct cats_CellsPool;
struct cats_WeightedLinksPool;

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
                        , const cats_HitData_t
                        );

/**\brief Returns number of points in layer
 *
 * Used mainly for debug purposes.
 *
 * \returns number of points added with `cats_layer_add_point()` for layer N
 * */
size_t cats_layer_n_points( const struct cats_Layers *, cats_LayerNo_t );

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

/**\brief Re-sets the pool for rentrant usage
 *
 * Useful for reentrant usage of the layers stack between switching the
 * connection topologies (change missing layers, apply doublet filtering or
 * lookup, etc).
 * */
void cats_cells_pool_reset( struct cats_Layers *
                          , struct cats_CellsPool *
                          , size_t softLimitCells
                          );

/**\brief Builds a pairwise connection graph to be evaluated on given layers
 *
 * This routine appends the given layers stack with connection information.
 * Connections are established totally, between neighbouring layers, without
 * doublet filtering or lookup, based only on the permission provided by
 * logically-discrete triplet filtering.
 *
 * \returns `CATSC_RC_EMPTY_GRAPH` if no connection can be performed for
 *          given layers
 * \returns `CATSC_ERROR_ALLOC_FAILURE_CELLS` at cell (link b/w hits)
 *          allocation failure.
 * \returns `CATSC_ERROR_ALLOC_FAILURE_BACKREF_COUNTER` on cell counter
 *          overflow (too many connections)
 * \returns `CATSC_ERROR_ALLOC_FAILURE_BACKREF` on back references list
 *          overflow (too many connections)
 * */
int cats_connect( struct cats_Layers * ls
                , struct cats_CellsPool * a
                , cats_Filter_t test_triplet
                , void * userData
                , cats_LayerNo_t nMissingLayers
                );

/**\brief Builds a pairwise connection graph to be evaluated on given layers
 *
 * This routine appends the given layers stack with connection information.
 * Connections are established totally, between neighbouring layers, without
 * doublet filtering or lookup, based only on the permission provided by
 * quasi-contiguous (floating point) triplet weighting.
 *
 * Return codes are the same as for `cats_connect()` plus following one:
 *
 * \returns `CATSC_ERROR_ALLOC_FAILURE_WEIGHTS` in case of (re)alocation
 *          failure for internal weights buffer.
 * */
int
cats_connect_w( struct cats_Layers * ls
              , struct cats_CellsPool * a
              , cats_WeightedFilter_t test_triplet
              , void * userData
              , cats_LayerNo_t nMissingLayers
              );

/**\brief Performs forward evaluation (1st stage) of CATS algorithm
 *
 * Constructs cellular automaton based on hits info and evaluates states till
 * convergence.
 *
 * Set `debugJSONStream` to disable debug dump for each iteration.
 *
 * \returns    0 if automaton is evaluated and ready to provide the tracks.
 * */
int cats_evaluate( struct cats_Layers *, FILE * debugJSONStream);

/**\brief Re-sets "visited" flags after previous `collect()` call
 *
 * Re-sets internal state flags used to mark visited cells. Should be used
 * between switching strategies.
 *
 * \note Only useful if connection topology can to be left intact between
 *       applying different strategies.
 * */
void reset_collection_flags( struct cats_CellsPool * );

/*                                            ________________________________
 * ________________________________________ / Unweighted collection strategies
 */


/**\brief Iterates over resulting connection graph visiting all the enumerated
 *        subsets, excessively (with sub-sequences).
 *
 * Will invoke given callback on each found track candidate, including
 * sub-sequences when they are permitted by non-weighted geometrical filter.
 *
 * Provides (efficient) equivalent to brute-forced connectionalgorithm,
 * where all combinations passing through the geometrical filter are
 * considered. Hit insertion order affects determinism: combinations
 * corresponging to hits inserted first will be considered first.
 *
 * This evaluation routine provides excessive candidates and can be used
 * for certain tracking strategies scrutinizing every possible track
 * combination permitted by some geometrical filter, or as a last resort for
 * cases of low efficiency.
 * */
int
cats_visit_dfs_excessive( struct cats_Layers * ls
                        , unsigned int minLength
                        , void (*callback)(const cats_HitData_t *, size_t, void *)
                        , void * userdata
                        );

/**\brief Iterates over resulting connection graph visiting all the enumerated
 *        subsets, permitting sub-sequences
 *
 * This is collecting routine will invoke given callback on each found track
 * candidate, including sub-sequences when they are permitted by non-weighted
 * geometrical filter. Hit insertion order
 * affects determinism: combinations corresponging to hits inserted first will
 * be considered first.
 *
 * This evaluation routine can be preferable for certain tracking scenarios
 * whith low multiplicity and low efficiency, albeit for high multiplicities
 * can create highly excessive redundancy.
 *
 * Comparing to excessive strategy it will not iterate over certain sequences
 * (TODO: elaborate the difference).
 * */
int
cats_visit_dfs_moderate( struct cats_Layers * ls
                       , unsigned int minLength
                       //, unsigned int nMissingLayers
                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                       , void * userdata
                       );

/**\brief Iterates over resulting connection graph visiting all the enumerated
 *        subsets, omitting sub-sequences
 *
 * This is collecting routine will invoke given callback on each found track
 * candidate, excluding their sub-sequences when they are permitted by
 * non-weighted geometrical filter. Hit insertion order
 * affects determinism: combinations corresponging to hits inserted first will
 * be considered first.
 *
 * This evaluation routine can be preferable for certain tracking scenarios
 * whith high multiplicity and low efficiency.
 * */
int
cats_visit_dfs_strict( struct cats_Layers * ls
                     , unsigned int minLength
                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                     , void * userdata
                     );

/**\brief Iterates over resulting connection graph visiting enumerated
 *        subsets, preferring longest sequences, by hit insertion order.
 *
 * This routine will traverse connection graph preferring longest sequences
 * over shorter ones (yet, still permitted by the unweighted geometrical
 * filter). For concurrent tracks, a hit inserted first will be preferred.
 *
 * This function might be useful for tracking strategies anticipating deviated
 * hits, like Kalman with deterministing annealing filter. It effectively reduces
 * combinatorics, yet still can pick up noisy hits and may incorrectly resolve
 * pile-up concurrency. Since it does not need the weighted filter, this case
 * can be still preferrable for tracking being sparse enough for some strict
 * yet efficient geometrical filter. Typical use case -- scarcely occupied
 * efficient but noisy detectors.
 *
 * \todo Has an issue of doubling endings; currently can be mitigated with
 *       `LongestUniqueTrackCollector`, but has to be fixed with introducing
 *       "visited" flag per hit
 * */
int
cats_visit_dfs_longest( struct cats_Layers * ls
                      , unsigned int minLength
                      , unsigned int nMissingLayers
                      , void (*callback)(const cats_HitData_t *, size_t, void *)
                      , void * userdata
                      );

/**\brief Iterates over resulting connection graph visiting enumerated
 *        subsets, preferring longest sequences, by hit insertion order.
 *
 * This routine will traverse connection graph preferring longest sequences
 * over shorter ones (yet, still permitted by the unweighted geometrical
 * filter). For concurrent tracks, a hit inserted first will be preferred.
 * Comparing to `cats_for_each_longest_track_candidate()` does not repeat
 * sub-sequences.
 *
 * Use case is similar to `cats_for_each_longest_track_candidate()`, but for
 * elaborated tracking procedures capable to combine tracklets from
 * reconstructed segments.
 * */
int
cats_visit_dfs_winning( struct cats_Layers * ls
                      , unsigned int minLength
                      , unsigned int nMissingLayers
                      , void (*callback)(const cats_HitData_t *, size_t, void *)
                      , void * userdata
                      );

/*                                              ______________________________
 * __________________________________________ / Weighted collection strategies
 */

int
cats_visit_dfs_excessive_w( struct cats_Layers * ls
                          , unsigned int minLength
                          , void (*callback)(const cats_HitData_t *, size_t, void *)
                          , void * userdata
                          );

/*
 * Strict weighted visiting strategy preferring closest hits over distant ones
 * even if distant ones are better in terms of weight.
 */

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  /* H_CATS_H */
