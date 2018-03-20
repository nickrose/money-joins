//===----------------------------------------------------------------------===//
//
// A regular open addressing table with linear probing. Inlined key/value pairs.
//
//===----------------------------------------------------------------------===//

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>              /* CPU_ZERO, CPU_SET */
#include <pthread.h>            /* pthread_* */
#include <string.h>             /* memset */
#include <stdio.h>              /* printf */
#include <stdlib.h>             /* memalign */
#include <sys/time.h>           /* gettimeofday */

#include "pmoney_join.h"
#include "pmj_params.h"         /* constant parameters */
#include "rdtsc.h"              /* startTimer, stopTimer */
#include "lock.h"               /* lock, unlock */
#include "cpu_mapping.h"        /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

#include "hash.h"
#include "generator.h"          /* numa_localize() */

#ifndef NEXT_POW_2
/** 
 *  compute the next number, greater than or equal to 32-bit unsigned v.
 *  taken from "bit twiddling hacks":
 *  http://graphics.stanford.edu/~seander/bithacks.html
 */
#define NEXT_POW_2(V)                           \
    do {                                        \
        V--;                                    \
        V |= V >> 1;                            \
        V |= V >> 2;                            \
        V |= V >> 4;                            \
        V |= V >> 8;                            \
        V |= V >> 16;                           \
        V++;                                    \
    } while(0)
#endif

// Debug msg logging method
#ifdef DEBUG
#define DEBUGMSG(COND, MSG, ...)                                    \
    if(COND) { fprintf(stderr, "[DEBUG] "MSG, ## __VA_ARGS__); }
#else
#define DEBUGMSG(COND, MSG, ...) 
#endif

#define POW2
#define PREFETCH

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

struct bucket_t {
  intkey_t hash;
  tuple_t  tuple;
};

typedef struct bucket_t bucket_t;

typedef struct hashtable {
  bucket_t  *buckets;
  uint32_t  num_buckets;
  uint32_t  hash_mask;
  uint32_t  skip_bits;
} hashtable_t;

static inline uint32_t alt_mod(uint32_t x, uint32_t n) {
  return ((uint64_t) x * (uint64_t) n) >> 32 ;
}

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
#ifdef POW2
  ht->num_buckets *= 2;
  NEXT_POW_2((ht->num_buckets));
#else
  ht->num_buckets = ((double)ht->num_buckets) * 1.2;
#endif

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE,
                     ht->num_buckets * sizeof(bucket_t))){
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  // Not an elegant way of passing whether we will numa-localize, but this
  // feature is experimental anyways
  /*
  if(numalocalize) {
    tuple_t * mem = (tuple_t *) ht->buckets;
    uint32_t ntuples = (ht->num_buckets*sizeof(vs_t))/sizeof(tuple_t);
    numa_localize(mem, ntuples, nthreads);
  }
  */

  fprintf(stderr, "HT: %.2lf KB, # Buckets: %u, LF: %.2lf\n",
          (ht->num_buckets*sizeof(bucket_t)) / 1024.0,
          ht->num_buckets,
          (double)nbuckets/ht->num_buckets);

  memset(ht->buckets, 0, ht->num_buckets * sizeof(bucket_t));
  ht->skip_bits = 0;
  ht->hash_mask = (ht->num_buckets - 1) << ht->skip_bits;
  *ppht = ht;
}

//===----------------------------------------------------------------------===//
// Destroy the given hash table
//===----------------------------------------------------------------------===//
static inline void destroy_hashtable(hashtable_t * ht) {
  free(ht->buckets);
  free(ht);
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//

static inline void build_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
#ifdef POW2
    uint32_t idx = hash & hashmask;
#else
    uint32_t idx = alt_mod(hash, ht->num_buckets);
#endif
    while (ht->buckets[idx].hash) {
#ifdef POW2
      idx = (idx + 1) & hashmask;
#else
      if (++idx == ht->num_buckets) idx = 0;
#endif
    }
    ht->buckets[idx].hash = hash;
    ht->buckets[idx].tuple = rel->tuples[i];
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  uint64_t matches = 0;
  //double avg_len = 0.0, min = 100000.0, max = 0.0;
#ifdef PREFETCH
  uint32_t pf_idx = 64;
#endif

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
#ifdef POW2
    uint32_t idx = hash & hashmask;
#else
    uint32_t idx = alt_mod(hash, ht->num_buckets);
#endif

    // Do prefetch if enabled
#ifdef PREFETCH
    if (pf_idx < rel->num_tuples) {
      uint32_t pf_hash = Hash(rel->tuples[pf_idx++].key);
#ifdef POW2
      uint32_t idx_prefetch = pf_hash & hashmask;
#else
      uint32_t idx_prefetch = alt_mod(pf_hash, ht->num_buckets);
#endif
      __builtin_prefetch(ht->buckets + idx_prefetch, 0, 1);
    }
#endif

    //uint32_t len = 1;
    while ((ht->buckets[idx].hash) &&
           (ht->buckets[idx].tuple.key != rel->tuples[i].key)) {
#ifdef POW2
      idx = (idx + 1) & hashmask;
#else
      if (++idx == ht->num_buckets) idx = 0;
#endif
      //len++;
    }
    //avg_len += len;
    //min = (len < min ? len : min);
    //max = (len > max ? len : max);
    if (ht->buckets[idx].hash) {
      matches++;
    }
  }
  //avg_len /= rel->num_tuples;
  //printf("Avg. len: %.2lf, min: %.2lf, max: %.2lf\n", avg_len, min, max);
  return matches;
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t PMJ_3(relation_t *relR, relation_t *relS, int nthreads) {

#ifndef NO_TIMING
  struct timeval start, end;
  uint64_t partition_time, build_time, probe_time;
#endif

  hashtable_t *ht;
  uint32_t nbuckets = relR->num_tuples;
  allocate_hashtable(&ht, nbuckets);

#ifdef PERF_COUNTERS
  PCM_initPerformanceMonitor(NULL, NULL);
  PCM_start();
#endif

#ifndef NO_TIMING
  gettimeofday(&start, NULL);
  startTimer(&build_time);
  startTimer(&probe_time); 
  partition_time = 0; // no partitioning
#endif

  //////////////////////////////////////
  // BUILD
  //////////////////////////////////////
  build_hashtable_st(ht, relR);

#if 1
  uint32_t occ = 0;
  for (uint32_t j = 0; j < ht->num_buckets; j++) {
    occ += (ht->buckets[j].hash != 0);
  }
  fprintf(stderr, "Actual load-factor: %.2lf \n",
          (double)occ/((double)ht->num_buckets));
#endif

#ifndef NO_TIMING
  stopTimer(&build_time);
#endif

#ifdef PERF_COUNTERS
  PCM_stop();
  PCM_log("========== Build phase profiling results ==========\n");
  PCM_printResults();
  PCM_start();
#endif

  //////////////////////////////////////
  // PROBE
  //////////////////////////////////////
  int64_t result = probe_hashtable_st(ht, relS);

#ifdef PERF_COUNTERS
    PCM_stop();
    PCM_log("========== Probe phase profiling results ==========\n");
    PCM_printResults();
    PCM_log("===================================================\n");
    PCM_cleanup();
#endif

#ifndef NO_TIMING
  stopTimer(&probe_time);
  gettimeofday(&end, NULL);
  // Now print the timing results
  print_timing(probe_time, build_time, partition_time,
               relR->num_tuples, relS->num_tuples, result, &start, &end);
#endif


  // Cleanup
  destroy_hashtable(ht);

  return result;
}
