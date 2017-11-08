//===----------------------------------------------------------------------===//
//
// Vectorized, bucket-chain hash table
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

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

typedef struct vs {
  tuple_t  tuple;
  uint32_t next;
} vs_t;

typedef struct hashtable {
  uint32_t  *buckets;
  vs_t      *values;
  uint32_t  first_empty;
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
  //ht->num_buckets <<= 1;

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE,
                     ht->num_buckets * sizeof(uint32_t))){
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  /* allocate hashtable valuespace cache line aligned */
  if (posix_memalign((void**)&ht->values, CACHE_LINE_SIZE,
                     (nbuckets+1) * sizeof(vs_t))){
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
  uint32_t next_sz = ht->num_buckets * sizeof(uint32_t);
  uint32_t vals_sz = (nbuckets+1) * sizeof(vs_t);
  fprintf(stderr, "HT: %.2lf KB (%u buckets, sizeof vs: %lu)\n",
          (next_sz+vals_sz) / 1024.0, ht->num_buckets, sizeof(vs_t));

  memset(ht->buckets, 0, ht->num_buckets * sizeof(uint32_t));
  memset(ht->values, 0, (nbuckets+1) * sizeof(vs_t));
  ht->hash_mask = (ht->num_buckets - 1);
  ht->first_empty = 1;
  *ppht = ht;
}

//===----------------------------------------------------------------------===//
// Destroy the given hash table
//===----------------------------------------------------------------------===//
static inline void destroy_hashtable(hashtable_t * ht) {
  free(ht->buckets);
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
    uint32_t pos = ht->first_empty;
    
    ht->values[pos].next = ht->buckets[idx];
    ht->buckets[idx] = pos;
    ht->values[pos].tuple = rel->tuples[i];
    ht->first_empty++;
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
#if 0
static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

  uint64_t matches = 0;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    uint32_t pos = ht->buckets[idx];
    if (pos) {
      do {
        if (rel->tuples[i].key == ht->values[pos].tuple.key) {
          matches++;
        }
      } while ( (pos = ht->values[pos].next) );
    }
  }
  return matches;
}
#endif

#if 0
static inline int64_t check_next(uint32_t n, tuple_t *tuples, hashtable_t *ht, uint32_t *pos, uint32_t *match) {
  uint32_t k = 0;
  for (uint32_t i = 0; i < n; i++) {
    if (tuples[match[i]].key != ht->values[pos[match[i]]].tuple.key) {
      pos[match[i]] = ht->values[pos[match[i]]].next; 
      match[k] = match[i];
      k += (pos[match[i]] != 0);
    }    
  }
  return k;
}

static inline int64_t vectorized_probe(hashtable_t *ht, tuple_t *tuples, uint32_t n, uint32_t *pos, uint32_t *match) {
  const uint32_t hashmask = ht->hash_mask;

  int64_t result = 0; 
  uint32_t k = 0;

  // Initial lookup
  for (uint32_t i = 0; i < n; i++) {
    uint32_t idx = Hash(tuples[i].key) & hashmask;
    //pos[i] = ht->buckets[idx]; 
    pos[k] = ht->buckets[idx]; 
    match[k] = i;
    k += (pos[i] != 0);
  }

#if 0
  for (uint32_t i = 0; i < n; i++) {
    pos[k] = pos[i];
    match[k] = i;
    k += (pos[i] != 0);
  }
#endif

  // Recheck
  while (k > 0) {
    k = check_next(k, tuples, ht, pos, match);
  }

  for (uint32_t i = 0; i < n; i++) {
    result += (pos[i] != 0);
  }

  return result;
}

static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  uint32_t pos[VECTOR_SIZE], matches[VECTOR_SIZE];
  uint64_t result = 0;

  for (uint32_t i = 0; i < rel->num_tuples; i += VECTOR_SIZE) {
    uint32_t sz = i + VECTOR_SIZE < rel->num_tuples ? VECTOR_SIZE : rel->num_tuples - i;
    result += vectorized_probe(ht, rel->tuples + i, sz, pos, matches);    
  }
  return result;
}
#endif

#if 1
static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  uint64_t matches = 0;
  uint32_t avg_len = 0, min = 100000, max = 0;

  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t idx = Hash(rel->tuples[i].key) & hashmask;
    uint32_t len = 0.0;
    for (uint32_t hit = ht->buckets[idx]; hit > 0; hit = ht->values[hit].next) {
      len++;
      matches += (rel->tuples[i].key == ht->values[hit].tuple.key);
    }
    avg_len += len;
    min = (len < min ? len : min);
    max = (len > max ? len : max);
  }
  printf("Avg. len: %.2lf, min: %u, max: %u\n",
         ((double)avg_len)/rel->num_tuples, min, max);
  return matches;
}
#endif

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
int64_t PMJ_2(relation_t *relR, relation_t *relS, int nthreads) {

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
    occ += (ht->buckets[j] != 0);
  }
  fprintf(stderr, "HT load-factor: %.2lf \n",
          (double)occ/((double)ht->num_buckets));
#endif

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
