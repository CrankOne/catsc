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
#define CATSC_ERROR_RUNTIME_LOGIC                   0x10  /* runtime error flag */
#define CATSC_RC_NO_SUCH_LAYER                      0x11  /* No layer of given number */
#define CATSC_RC_EMPTY_GRAPH                        0x12  /* Graph has no connections (not actually erorr) */
#define CATSC_RC_BAD_LAYERS_NO                      0x13  /* Layers number is too small */

#ifndef CATS_HIT_DATA_TYPE
#   define CATS_HIT_DATA_TYPE const void*
#endif

/** Type of hit data payload */
typedef CATS_HIT_DATA_TYPE cats_HitData_t;

/** Filter callback type */
typedef int (*cats_Filter_t)( cats_HitData_t
                            , cats_HitData_t
                            , cats_HitData_t
                            , void * );

/** Weighted filter callback type */
typedef double (*cats_WeightedFilter_t)( cats_HitData_t
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

/**\brief Re-sets the pool for rentrant usage*/
void cats_cells_pool_reset( struct cats_Layers *
                          , struct cats_CellsPool *
                          , size_t softLimitCells
                          );

/**\brief Performs forward evaluation (1st stage) of CATS algorithm
 *
 * Constructs cellular automaton based on hits info and evaluates states till
 * convergence.
 *
 * Set `debugJSONStream` to disable debug dump for each iteration.
 *
 * \returns    0 if automaton is evaluated and ready to provide the tracks.
 * \returns -101 at cell (link b/w hits) allocation failure.
 * \returns -102 at cell neighbor reference allocation failure.
 * */
int cats_evaluate( struct cats_Layers *
                 , struct cats_CellsPool *
                 , cats_Filter_t test_triplet
                 , void * userData
                 , unsigned int nMissingLayers
                 , FILE * debugJSONStream
                 );

/**\brief Re-sets "visited" flags after previous `collect()` call
 *
 * Re-sets internal state flags used to mark visited cells.
 * */
void reset_collection_flags( struct cats_CellsPool * );

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
void
cats_for_each_track_candidate_excessive( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdata
                                       );

/**\brief Iterates over resulting connection graph visiting all the enumerated
 *        subsets, permitting duplicating sub-sequences due to inefficiences
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
 * */
void
cats_for_each_track_candidate( struct cats_Layers * ls
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
 *
 * This function might be useful for tracking strategies anticipating deviated
 * hits, like Kalman with deterministing annealing filter. It effectively reduces
 * combinatorics, yet still can pick up nosiy hits and may incorrectly resolve
 * pile-up concurrency. Since it does not need the weighted filter, this case
 * can be still preferrable for tracking being sparse enough for some strict
 * yet efficient geometrical filter.
 *
 * \todo Has an issue of doubling endings; currently can be mitigated with
 *       `LongestUniqueTrackCollector`, but has to be fixed with introducing
 *       "visited" flag per hit
 * */
void
cats_for_each_longest_track_candidate( struct cats_Layers * ls
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
void
cats_for_each_winning_track_candidate( struct cats_Layers * ls
                                     , unsigned int minLength
                                     , unsigned int nMissingLayers
                                     , void (*callback)(const cats_HitData_t *, size_t, void *)
                                     , void * userdata
                                     );

/**\brief Iterates over resulting graph, picking up weighted connections
 *
 * Choice among hits with concurrent weights will be resolved in favour of the
 * ones (with highest weight), according to given weighting filter. If two
 * weights are the same, first inserted of closest layer is preferred. Hit
 * insertion order should not affect determinism.
 *
 * Effectively reduces combinatorics, yet loses some candidates with
 * systematically loose weight (e.g. misaligned ones). This collecting routine
 * can be preferred for systems with high detector redundancy, on high
 * occupancies.
 * */
void
cats_for_each_track_candidate_w( struct cats_Layers * ls
                               , unsigned int minLength
                               , unsigned int nMissingLayers
                               , cats_WeightedFilter_t
                               , void * userdataFilter
                               , void (*callback)(const cats_HitData_t *, size_t, void *)
                               , void * userdataCollect
                               );

/**\brief Iterates over resulting graph, picking up longest weighted connections
 *
 * Concurrent cells will be resolved in favor of the best among closest ones,
 * according to weighting function. If two weights are the same, first inserted
 * is preferred. Hit insertion order should not affect determinism.
 *
 * Effectively reduces combinatorics, yet loses some candidates with
 * systematically loose weight (e.g. misaligned ones). This collecting routine
 * is preferred for systems with low redundancy on moderate occuancies.
 * */
void
cats_for_each_longest_track_candidate_w( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers
                                       , cats_WeightedFilter_t
                                       , void * userdataFilter
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdataCollect
                                       );

/**\brief Iterates over resulting graph, picking up longest weighted
 *        connections, never passing through the same segments twice
 *
 * Concurrent cells will be resolved in favor of the best among closest ones,
 * according to weighting function. If two weights are the same, first inserted
 * is preferred. Hit insertion order should not affect determinism. This
 * implementation follows `winner takes all` rule, producing candidates without
 * shared segments which can hide pile-up events, yet providing best possible
 * combinatorial reduction.
 *
 * Effectively reduces combinatorics, yet loses some candidates with
 * systematically loose weight (e.g. misaligned ones). This collecting routine
 * is preferred for systems with low redundancy on moderate occuancies.
 * */
void
cats_for_each_winning_track_candidate_w( struct cats_Layers * ls
                                       , unsigned int minLength
                                       , unsigned int nMissingLayers
                                       , cats_WeightedFilter_t
                                       , void * userdataFilter
                                       , void (*callback)(const cats_HitData_t *, size_t, void *)
                                       , void * userdataCollect
                                       );

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  /* H_CATS_H */
