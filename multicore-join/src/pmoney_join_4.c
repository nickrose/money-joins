//===----------------------------------------------------------------------===//
//
// Open addressing, linear probing using bitmap status
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

#include "pmoney_join_4.h"
#include "pmj_params.h"         /* constant parameters */
#include "rdtsc.h"              /* startTimer, stopTimer */
#include "lock.h"               /* lock, unlock */
#include "cpu_mapping.h"        /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

#include "hash.h"
#include "generator.h"          /* numa_localize() */

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

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

#ifndef HASH
#define HASH(X, MASK, SKIP) (((Hash(X)) & MASK) >> SKIP)
#endif

// Debug msg logging method
#ifdef DEBUG
#define DEBUGMSG(COND, MSG, ...)                                    \
    if(COND) { fprintf(stderr, "[DEBUG] "MSG, ## __VA_ARGS__); }
#else
#define DEBUGMSG(COND, MSG, ...) 
#endif

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

typedef struct hashtable {
  uint8_t   *flags;
  tuple_t   *values;
  uint32_t  num_buckets;
  uint32_t  hash_mask;
  uint32_t  skip_bits;
} hashtable_t;

#define SET_BIT(bitmap, idx) bitmap[idx >> 3] |= 1 << (idx & 7)
#define IS_SET(bitmap, idx) (bitmap[idx >> 3] & (1 << (idx & 7)))

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->flags, CACHE_LINE_SIZE,
                     ht->num_buckets / 8 * sizeof(uint8_t))){
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }
  if (posix_memalign((void**)&ht->values, CACHE_LINE_SIZE,
                     ht->num_buckets * sizeof(tuple_t))){
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

  memset(ht->flags, 0, ht->num_buckets / 8 * sizeof(uint8_t));
  memset(ht->values, 0, ht->num_buckets * sizeof(tuple_t));
  ht->skip_bits = 0;
  ht->hash_mask = (ht->num_buckets - 1) << ht->skip_bits;
  *ppht = ht;
}

//===----------------------------------------------------------------------===//
// Destroy the given hash table
//===----------------------------------------------------------------------===//
static inline void destroy_hashtable(hashtable_t * ht) {
  free(ht->flags);
  free(ht->values);
  free(ht);
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//

static inline void build_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    uint32_t step = 0;
    do {
      if (IS_SET(ht->flags, idx) == 0) {
        ht->values[idx] = rel->tuples[i]; 
        SET_BIT(ht->flags, idx);
        break;
      }
      idx = (idx + 1) & hashmask;
    } while (step++ < ht->num_buckets);
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
#if 0
static inline int64_t do_probe_hashtable_st(hashtable_t *ht, tuple_t *tuples, uint32_t sz) {
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

  uint32_t matches[VECTOR_SIZE];
  uint32_t pivot[VECTOR_SIZE];
  uint32_t p = 0, m = 0;

  for (uint32_t i = 0; i < sz; i++) {
    uint32_t idx = HASH(tuples[i].key, hashmask, skipbits);
    matches[p] = idx;
    pivot[p] = i;
    p += (IS_SET(ht->flags, idx) != 0);
  }
  for (uint32_t i = 0; i < p; i++) {
    if (tuples[pivot[i]].key == ht->values[matches[i]].key) {
      m++;
    }
  }
  return m;

}

static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {

  uint64_t matches = 0;

  for (uint32_t i = 0; i < rel->num_tuples; i += VECTOR_SIZE) {
    uint32_t sz = i + VECTOR_SIZE < rel->num_tuples ? VECTOR_SIZE : rel->num_tuples - i;
    matches += do_probe_hashtable_st(ht, rel->tuples + i, sz);
  }
  return matches;
}
#endif

static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

  uint64_t matches = 0;
  uint32_t pf_idx = 16;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
#if 1
    if (pf_idx < rel->num_tuples) {
      uint32_t idx_pf = HASH(rel->tuples[pf_idx++].key, hashmask, skipbits);
      __builtin_prefetch(ht->values + idx_pf, 0, 1);
    }
#endif
    uint32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    if (IS_SET(ht->flags, idx)) {
      while (IS_SET(ht->flags, idx)) {
        if (ht->values[idx].key == rel->tuples[i].key) {
          matches++;
          break;
        }
        idx = (idx + 1) & hashmask;
      }
    }
  }
  return matches;
}

//===----------------------------------------------------------------------===//
// Print out the execution time statistics of the join
//===----------------------------------------------------------------------===//
static inline void print_timing(uint64_t total, uint64_t build, uint64_t part,
                                uint64_t num_build, uint64_t num_probe, int64_t result,
                                struct timeval *start, struct timeval *end) {
  // General
  uint64_t num_tuples = num_probe + num_build;
  double diff_msec = (((*end).tv_sec*1000L + (*end).tv_usec/1000L)
                      - ((*start).tv_sec*1000L+(*start).tv_usec/1000L));
  double cyclestuple = (double)total / (double)(num_tuples);

  // Throughput in million-tuples-per-sec
  double throughput = (double)num_tuples / (double)diff_msec / 1000.0;

  // Probe info
  uint64_t probe_cycles = total - build;
  double probe_cpt = (double)probe_cycles / (double)num_probe;
  double probe_usec = ((double)probe_cycles / (double)total) * diff_msec;

  // Build info
  uint64_t build_cycles = build - part;
  double build_cpt = (double)build_cycles / (double)num_probe;
  double build_usec = ((double)build_cycles / (double)total) * diff_msec;

  // Part
  double part_cpt = (double)part / (double)num_probe;
  double part_usec = ((double)part / (double)total) * diff_msec;

  fprintf(stderr, 
          "RESULTS: %lu, Runtime: %.2lf ms, Throughput: %.2lf mtps, " 
          "Probe: %.2lf ms (%.2lf CPT), Build: %.2lf ms (%.2lf CPT), "
          "Part: %.2lf ms (%.2lf CPT), CPT: %.4lf\n",
          result, diff_msec, throughput, probe_usec, 
          probe_cpt, build_usec, build_cpt, part_usec, part_cpt, cyclestuple);
  fprintf(stderr, "Total cycles: %lu\n", total);
  fflush(stderr);
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t PMJ_4(relation_t *relR, relation_t *relS, int nthreads) {

#ifndef NO_TIMING
  struct timeval start, end;
  uint64_t partition_time, build_time, probe_time;
#endif

  hashtable_t *ht;
  uint32_t nbuckets = 2 * relR->num_tuples;
  size_t mem_before = get_memory_usage_bytes();
  allocate_hashtable(&ht, nbuckets);
  size_t mem_after = get_memory_usage_bytes();
  fprintf(stderr, "HT: %.2lf MB\n",
          ((double)mem_after-mem_before)/1024.0/1024.0);

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

#ifdef DEBUG
  uint32_t occ = 0;
  for (uint32_t j = 0; j < ht->num_buckets; j++) {
    occ += (ht->buckets[j] != 0);
  }
  fprintf(stderr, "HT load-factor: %.2lf \%\n", (double)occ/(double)ht->num_buckets);
  fprintf(stderr, "Num allocs: %llu\n", num_allocs);
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
