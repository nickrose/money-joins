
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

#define BUCKET_SIZE 8

typedef struct bucket {
  uint32_t hashes[BUCKET_SIZE];
  intkey_t keys[BUCKET_SIZE];
  value_t vals[BUCKET_SIZE];
  struct bucket *next;
} __attribute__((__packed__))  bucket_t;

typedef struct hashtable {
  bucket_t      *buckets;
  bucket_t      *pool;
  bucket_t      *pool_pos;
  uint32_t      num_buckets;
  uint32_t      hash_mask;
} hashtable_t;

// Fast alternative to modulo from Daniel Lemire
uint32_t alt_mod(uint32_t x, uint32_t n) {
  return ((uint64_t) x * (uint64_t) n) >> 32 ;
}

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nelems) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = ((double)nelems)*1.2;
  //ht->num_buckets = nelems;
  //NEXT_POW_2((ht->num_buckets));
  ht->num_buckets /= BUCKET_SIZE;
 
  uint64_t required_bytes = ht->num_buckets * sizeof(bucket_t);
  if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE, required_bytes)) {
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  uint64_t overflow = (nelems / BUCKET_SIZE) * sizeof(bucket_t);
  if (posix_memalign((void**)&ht->pool, CACHE_LINE_SIZE, overflow)) {
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

  fprintf(stderr, "HT: %.2lf KB (%u buckets of size %u each)\n",
          required_bytes/1024.0, ht->num_buckets, sizeof(bucket_t));

  memset(ht->buckets, 0, required_bytes);
  memset(ht->pool, 0, overflow);
  ht->pool_pos = ht->pool;
  ht->hash_mask = ht->num_buckets - 1;
  *ppht = ht;
}

//===----------------------------------------------------------------------===//
// Destroy the given hash table
//===----------------------------------------------------------------------===//
static inline void destroy_hashtable(hashtable_t * ht) {
  free(ht->buckets);
  free(ht->pool);
  free(ht);
}

void dump_bucket(bucket_t *b) {
  uint32_t nelems = 0;
  for (uint32_t i = 0; i < BUCKET_SIZE; i++) {
    if (i != 0) fprintf(stderr, ",");
    fprintf(stderr, "hashes[%u]=%u", i, b->hashes[i]);
    if (b->hashes[i] != 0) nelems++;
  }
  fprintf(stderr, ", #elems = %u\n", nelems);
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//

static inline void build_hashtable_st(hashtable_t *ht, relation_t *rel) {
  const uint32_t hashmask = ht->hash_mask;

  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
    bucket_t *b = ht->buckets + /*(hash & hashmask)*/alt_mod(hash, ht->num_buckets);
    bool found = false;
    while (b && !found) {
      for (uint32_t j = 0; j < BUCKET_SIZE; j++) {
        if (b->hashes[j] == 0) {
          b->hashes[j] = hash;
          b->keys[j] = rel->tuples[i].key;
          b->vals[j] = rel->tuples[i].payload;
          found = true;
          break;
        }
      }
      b = b->next;
    }
    if (!found) {
      //fprintf(stderr, "allocating overflow\n");
      bucket_t *newb = ht->pool_pos;
      bucket_t *b = ht->buckets + (hash & hashmask);
      bucket_t *curr_next = b->next;

      newb->next = curr_next;
      b->next = newb;

      ht->pool_pos++;
      newb->hashes[0] = hash;
      newb->keys[0] = rel->tuples[i].key;
      newb->vals[0] = rel->tuples[i].payload;
    }
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable_st(hashtable_t *ht, relation_t *rel, relation_t *rel_build) {
  const uint32_t hashmask = ht->hash_mask;

  uint64_t matches = 0;
  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    uint32_t hash = Hash(rel->tuples[i].key);
    bucket_t *b = ht->buckets + /*(hash & hashmask)*/alt_mod(hash, ht->num_buckets);
    for (; b; b = b->next) {
      uint32_t j = 0;
      for (; j < BUCKET_SIZE; j++) {
        if (b->hashes[j] == hash && b->keys[j] == rel->tuples[i].key) break;
      }
      if (j < BUCKET_SIZE) {
        matches++;
        break;
      }
    }
  }
  return matches;
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t PMJ_9(relation_t *relR, relation_t *relS, int nthreads) {

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
    for (uint32_t x = 0; x < BUCKET_SIZE; x++) {
      occ += (ht->buckets[j].hashes[x] != 0);
    }
  }
  fprintf(stderr, "HT load-factor: %.2lf \n",
          (double)occ/((double)ht->num_buckets/BUCKET_SIZE));
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
  int64_t result = probe_hashtable_st(ht, relS, relR);

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
