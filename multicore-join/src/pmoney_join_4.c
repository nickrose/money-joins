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

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

typedef struct hashtable {
  uint8_t   *flags;
  tuple_t   *values;
  uint32_t  num_buckets;
  uint32_t  hash_mask;
} hashtable_t;

#define POW2

#define SET_BIT(bitmap, idx) bitmap[idx >> 3] |= 1 << (idx & 7)
//#define IS_SET(bitmap, idx) ((bitmap[idx >> 3] & (1 << (idx & 7))) != 0)
#define IS_SET(bitmap, idx) ((bitmap[idx >> 3] & (1 << (idx & 7))))

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));
  //ht->num_buckets <<= 1;

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->flags, CACHE_LINE_SIZE,
                     (ht->num_buckets / 8) * sizeof(uint8_t))){
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
  uint32_t flag_sz = (ht->num_buckets / 8) * sizeof(uint8_t);
  uint32_t buc_sz = ht->num_buckets * sizeof(tuple_t);
  fprintf(stderr, "HT: %.2lf KB (bitmap: %.2lf KB, values: %.2lf KB), # Buckets: %u, LF: %.2lf\n",
          (flag_sz+buc_sz)/1024.0, flag_sz / 1024.0, buc_sz / 1024.0,
          ht->num_buckets,
          (double)nbuckets/ht->num_buckets);

  memset(ht->flags, 0, ht->num_buckets / 8 * sizeof(uint8_t));
  memset(ht->values, 0, ht->num_buckets * sizeof(tuple_t));
  ht->hash_mask = (ht->num_buckets - 1);
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

  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = Hash(rel->tuples[i].key) & hashmask;
    while (IS_SET(ht->flags, idx)) {
      idx = (idx + 1) & hashmask;
    }
    SET_BIT(ht->flags, idx);
    ht->values[idx] = rel->tuples[i];
  }
}

static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  uint64_t matches = 0;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
    uint32_t idx = hash & hashmask;
#if 0
    // Fast existence check, not actually probe since we don't check for
    // key or hash equality ....
    matches += (IS_SET(ht->flags, idx));
#else
    while (IS_SET(ht->flags, idx) &&
           ht->values[idx].key != rel->tuples[i].key) {
      idx = (idx + 1) & hashmask;
    }
    if (IS_SET(ht->flags, idx)) {
      matches++;
    }
#endif
  }
  return matches;
}

//===----------------------------------------------------------------------===//
// Print out the execution time statistics of the join
//===----------------------------------------------------------------------===//
static inline void print_timing(uint64_t total, double probe_clock, uint64_t build, uint64_t part,
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
  double probe_ns = (double)probe_clock/ CLOCKS_PER_SEC / (double)num_probe * 1e9;

  // Build info
  uint64_t build_cycles = build - part;
  double build_cpt = (double)build_cycles / (double)num_probe;
  double build_usec = ((double)build_cycles / (double)total) * diff_msec;

  // Part
  double part_cpt = (double)part / (double)num_probe;
  double part_usec = ((double)part / (double)total) * diff_msec;

  fprintf(stderr, 
          "RESULTS: %lu, Runtime: %.2lf ms, Throughput: %.2lf mtps, " 
          "Probe: %.2lf ms (%.2lf CPT, %.2lf ns/probe), Build: %.2lf ms (%.2lf CPT), "
          "Part: %.2lf ms (%.2lf CPT), CPT: %.4lf\n",
          result, diff_msec, throughput,
          probe_usec, probe_cpt, probe_ns,
          build_usec, build_cpt,
          part_usec, part_cpt,
          cyclestuple);
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
  uint32_t nbuckets = relR->num_tuples;
  size_t mem_before = get_memory_usage_bytes();
  allocate_hashtable(&ht, nbuckets);
  size_t mem_after = get_memory_usage_bytes();

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
  auto start_clock = clock();
  int64_t result = probe_hashtable_st(ht, relS);
  auto end_clock = clock();

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
  print_timing(probe_time, (double)(end_clock-start_clock), build_time, partition_time,
               relR->num_tuples, relS->num_tuples, result, &start, &end);
#endif


  // Cleanup
  destroy_hashtable(ht);

  return result;
}
