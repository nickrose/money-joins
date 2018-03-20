//===----------------------------------------------------------------------===//
//
// HyPER style joins
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
#include <immintrin.h>

#include "hyper_join.h"
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

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

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

typedef struct entry {
  tuple_t *tuple;
  struct entry *next; 
} entry_t;

typedef struct hashtable {
  entry_t   **buckets;

  char      *entry_mem_pool;
  char      *curr_mem_pool_pos;

  uint32_t  num_buckets;
  uint32_t  hash_mask;
} hashtable_t;

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));

  // Total number of bytes needed to store all entries
  uint64_t directory_size = ht->num_buckets * sizeof(entry_t **);
  uint64_t entry_pool_size = nbuckets * sizeof(entry_t);

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE, directory_size)) {
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  if (posix_memalign((void**)&ht->entry_mem_pool, CACHE_LINE_SIZE, entry_pool_size)){
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

  fprintf(stderr, "HT: %.2lf KB (%u buckets, %.2lf KB directory, %.2lf entry pool), LF: %.2lf\n",
          (directory_size + entry_pool_size) / 1024.0, ht->num_buckets,
          (double)directory_size/1024.0, (double)entry_pool_size/1024.0,
          (double) nbuckets / (double) ht->num_buckets / 1.0);

  memset(ht->buckets, 0, ht->num_buckets * sizeof(entry_t **));
  memset(ht->entry_mem_pool, 0, entry_pool_size);

  ht->curr_mem_pool_pos = ht->entry_mem_pool;
  ht->hash_mask = ht->num_buckets - 1;
  *ppht = ht;
}

//===----------------------------------------------------------------------===//
// Destroy the given hash table
//===----------------------------------------------------------------------===//
static inline void destroy_hashtable(hashtable_t * ht) {
  free(ht->buckets);
  free(ht->entry_mem_pool);
  free(ht);
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//

static inline void build_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = Hash(rel->tuples[i].key) & hashmask;
    entry_t *e = (entry_t *)ht->curr_mem_pool_pos;
    e->tuple = rel->tuples + i;
    e->next = ht->buckets[idx];
    ht->buckets[idx] = e;
    ht->curr_mem_pool_pos += sizeof(entry_t);
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  uint64_t matches = 0/*, steps = 0, curr_min = 1000000, curr_max = 0*/;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = Hash(rel->tuples[i].key) & hashmask;
    entry_t *e = ht->buckets[idx];
    //uint32_t step = 0;
    if (e) {
      do {
        if (e->tuple->key == rel->tuples[i].key) {
          matches++;
          break;
        }
        e = e->next;
        //step++; 
      } while (e);
    }
    //curr_max = max(curr_max, step);
    //curr_min = min(curr_min, step);
    //steps += step;
  }
  //double avg_step = (double)steps/rel->num_tuples;
  //fprintf(stderr, "Avg step: %.1lf, min: %lu, max: %lu\n",
  //        avg_step, curr_min, curr_max);
  return matches;
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t Hyper(relation_t *relR, relation_t *relS, int nthreads) {

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
