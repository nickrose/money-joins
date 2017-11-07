
// A regular linear probing HT

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
  uint32_t  mpl;
  uint32_t  hash_mask;
  uint32_t  skip_bits;
} hashtable_t;

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));
  //ht->num_buckets = ht->num_buckets << 1;

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

  memset(ht->buckets, 0, ht->num_buckets * sizeof(bucket_t));
  ht->mpl = 0;
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
    uint32_t idx = hash & hashmask;
    //uint32_t step = 1;
    while (ht->buckets[idx].hash) {
      idx = (idx + 1) & hashmask;
      //step <<= 1;
      //idx = (idx + step) & hashmask;
    }
    //printf("Key %lu -> idx %u\n", rel->tuples[i].key, idx);
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
  //float len = 0.0, mi = 10000000.0, ma = 0.0;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
    uint32_t idx = hash & hashmask;
    //float ilen = 1.0;
    //uint32_t step = 1;
    while ((ht->buckets[idx].hash) &&
           (ht->buckets[idx].tuple.key != rel->tuples[i].key)) {
      idx = (idx + 1) & hashmask;
      //ilen++;
    }
    //len += ilen;
    //mi = ilen < mi ? ilen : mi;
    //ma = ilen > ma ? ilen : ma;
    if (ht->buckets[idx].hash) {
      matches++;
    }
  }
  //len /= matches;
  //fprintf(stderr, "avg: %.2lf, min: %.2lf, max: %.2lf\n", len, mi, ma);
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
  fflush(stderr);
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

#if 1
  uint32_t occ = 0;
  for (uint32_t j = 0; j < ht->num_buckets; j++) {
    occ += ht->buckets[j].hash != 0;
  }
  fprintf(stderr, "HT load-factor: %.2lf \n",
          (double)occ/((double)ht->num_buckets*BUCKET_SIZE));
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
