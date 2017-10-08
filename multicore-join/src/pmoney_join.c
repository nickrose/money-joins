
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
//#include "npj_types.h"          /* bucket_t, hashtable_t, bucket_buffer_t */
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
#define HASH(X, MASK, SKIP) (((X) & MASK) >> SKIP)
#endif

// Debug msg logging method
#ifdef DEBUG
#define DEBUGMSG(COND, MSG, ...)                                    \
    if(COND) { fprintf(stderr, "[DEBUG] "MSG, ## __VA_ARGS__); }
#else
#define DEBUGMSG(COND, MSG, ...) 
#endif

//static uint64_t num_allocs = 0;

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

// Added the linear probing hash table when Instruments/Vtune/Perf said the
// probe loop was too hot ...
struct bucket_t {
  uint32_t flag;
  tuple_t tuple;
} /*__attribute__ ((aligned(CACHE_LINE_SIZE)))*/;
typedef struct bucket_t bucket_t;

struct hashtable_t {
  bucket_t *buckets;
  int32_t   num_buckets;
  uint32_t  hash_mask;
  uint32_t  skip_bits;
};
typedef struct hashtable_t hashtable_t;

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline void allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets) {
  hashtable_t * ht = (hashtable_t*) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE,
                     ht->num_buckets * sizeof(bucket_t))){
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  // Not an elegant way of passing whether we will numa-localize, but this
  // feature is experimental anyway.s
  if(numalocalize) {
    tuple_t * mem = (tuple_t *) ht->buckets;
    uint32_t ntuples = (ht->num_buckets*sizeof(bucket_t))/sizeof(tuple_t);
    numa_localize(mem, ntuples, nthreads);
  }

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

struct queue_t {
  uint64_t *positions;
};
typedef struct queue_t queue_t;

// Per-thread structure
struct arg_t {
  uint32_t  nthreads;  // The number of threads running the join
  uint32_t  tid;       // The tid of this thread

  uint32_t fanout;
  hashtable_t **hts;   // All the hashtables

  relation_t *build;
  relation_t *probe;
};
typedef struct arg_t arg_t;

struct hashes_t {
  int32_t hash;
  int32_t pos;
};
typedef struct hashes_t hashes_t;

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//

static inline void build_hashtable_st(arg_t *arg) {
  uint32_t hist[arg->fanout];

  hashes_t hashes[VECTOR_SIZE];
  hashes_t hashes_out[VECTOR_SIZE];
  memset(hashes, 0, sizeof(hashes));
  memset(hashes_out, 0, sizeof(hashes_out));

  uint32_t fanout = arg->fanout;
  uint32_t fanout_mask = fanout - 1;

  for (uint32_t i = 0; i < arg->build->num_tuples; i += VECTOR_SIZE) {
    uint32_t start = i;
    uint32_t end = arg->build->num_tuples < start + VECTOR_SIZE 
                     ? arg->build->num_tuples - start : VECTOR_SIZE;
    
    // First partition the input
    memset(hist, 0, sizeof(hist));

    // Hist
    for (uint32_t k = 0; k < end; k++) {
      hashes[k].hash = Hash(arg->build->tuples[start+k].key);
      hashes[k].pos = start + k;
      hist[hashes[k].hash & fanout_mask]++;
    }

    // Prefix sum
    for (uint32_t k = 0, sum = 0; k < fanout; k++) {
      int32_t tmp = hist[k];
      hist[k] = sum;
      sum += tmp;
    }
    hist[0] = 0;

    // Re-arrange
    for (uint32_t k = 0; k < end; k++) {
      uint32_t dest_idx = hist[hashes[k].hash & fanout_mask];
      hashes_out[dest_idx] = hashes[k];
      hist[hashes[k].hash & fanout_mask]++;
    }

    // Build
    for (uint32_t k = 0; k < end; k++) {
      hashtable_t *ht = arg->hts[hashes_out[k].hash & fanout_mask];
      uint32_t idx = (hashes_out[k].hash & ht->hash_mask) >> ht->skip_bits;
      
      uint32_t step = 0;
      while (1) {
        bucket_t *curr = &ht->buckets[idx];
        if (curr->flag == 0) {
          curr->tuple = arg->build->tuples[hashes_out[k].pos];
          curr->flag = 1;
          break;
        }
        idx = (idx + 1) & (ht->num_buckets - 1);
        step++;
      }
      if (step >= ht->num_buckets) {
        fprintf(stderr, "STEP: %u\n", step);
        exit(1);
      }
    }
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable_st(arg_t *arg) {
  uint32_t hist[arg->fanout];

  hashes_t hashes[VECTOR_SIZE];
  hashes_t hashes_out[VECTOR_SIZE];
  memset(hashes, 0, sizeof(hashes));
  memset(hashes_out, 0, sizeof(hashes_out));

  uint32_t fanout = arg->fanout;
  uint32_t fanout_mask = fanout - 1;

  uint64_t matches = 0;
  for (uint32_t i = 0; i < arg->probe->num_tuples; i += VECTOR_SIZE) {
    uint32_t start = i;
    uint32_t end = arg->probe->num_tuples < start + VECTOR_SIZE 
                     ? arg->probe->num_tuples - start : VECTOR_SIZE;
    
    // First partition the input
    memset(hist, 0, sizeof(hist));

    // Hist
    for (uint32_t k = 0; k < end; k++) {
      hashes[k].hash = Hash(arg->probe->tuples[start+k].key);
      hashes[k].pos = start + k;
      hist[hashes[k].hash & fanout_mask]++;
    }
    /*
    for (uint32_t k = 0; k < fanout; k++) {
      if (k % 10 == 0) fprintf(stderr, "\n");
      else if (k != 0) fprintf(stderr, ", ");
      fprintf(stderr, "%u -> %u", k, hist[k]);
    }
    fprintf(stderr, "\n");
    */

    // Prefix sum
    for (uint32_t k = 0, sum = 0; k < fanout; k++) {
      int32_t tmp = hist[k];
      hist[k] = sum;
      sum += tmp;
    }
    hist[0] = 0;

    // Re-arrange
    for (uint32_t k = 0; k < end; k++) {
      uint32_t dest_idx = hist[hashes[k].hash & fanout_mask];
      hashes_out[dest_idx] = hashes[k];
      hist[hashes[k].hash & fanout_mask]++;
    }

    // Probe
    for (uint32_t k = 0; k < end; k++) {
      hashtable_t *ht = arg->hts[hashes_out[k].hash & fanout_mask];
      uint32_t idx = HASH(hashes_out[k].hash, ht->hash_mask, ht->skip_bits);

      bucket_t *b = ht->buckets + idx;

      tuple_t *probe_tuple = &arg->probe->tuples[hashes_out[k].pos];

      //uint32_t m = matches;
      while (b->flag == 1) {
        if (b->tuple.key == probe_tuple->key) {
          matches++;
          break;
        }
        idx = (idx + 1) & (ht->num_buckets-1);
        b = ht->buckets + idx;
      }
      //fprintf(stderr, "%u has %u matches\n", probe_tuple->key, matches-m);
#if 0
      // Loop over all buckets till we find the end of the linked list
      do {
        for (uint32_t j = 0; j < b->count; j++) {
          if (probe_tuple->key == b->tuples[j].key) {
            matches++; 
          }
        }

        // Follow overflow pointer
        b = b->next;
      } while (b);
#endif
    }
  }
  return matches;
}

//===----------------------------------------------------------------------===//
// Print out the execution time statistics of the join
//===----------------------------------------------------------------------===//
static inline void print_timing(uint64_t total, uint64_t build, uint64_t part,
                                uint64_t numtuples, int64_t result,
                                struct timeval *start, struct timeval *end)
{
    double diff_msec = (((*end).tv_sec*1000L + (*end).tv_usec/1000L)
                        - ((*start).tv_sec*1000L+(*start).tv_usec/1000L));
    double cyclestuple = total;
    cyclestuple /= numtuples;
    uint64_t probe_cycles = total - build;
    double probe_usec = ((double)probe_cycles / (double)total) * diff_msec;
    uint64_t build_cycles = build - part;
    double build_usec = ((double)build_cycles / (double)total) * diff_msec;
    double part_usec = ((double)part / (double)total) * diff_msec;

    fprintf(stderr, 
            "RESULTS -> Runtime: %.2lf (%llu cycles), Probe: %.2lf (%llu cycles), "
            "Build: %.2lf (%llu cycles), Part: %.2lf (%llu cycles), CPT: %.4lf\n",
            diff_msec, total, probe_usec, probe_cycles, build_usec, build_cycles, 
            part_usec, part, cyclestuple);
    fflush(stderr);
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t PMJ(relation_t *relR, relation_t *relS, int nthreads) {
#ifndef NO_TIMING
  struct timeval start, end;
  uint64_t partition_time, build_time, probe_time;
#endif

  const uint32_t cache_size = 1 << 24; // 16MB cache
  const uint32_t num_tuples_in_cache = cache_size / sizeof(bucket_t);
  uint32_t fanout = relR->num_tuples / num_tuples_in_cache;
  NEXT_POW_2((fanout));
  //fanout = fanout < 2 ? 2 : fanout;
  uint32_t FANOUT_BITS = 64-__builtin_clzl(fanout-1);
  const uint32_t nbuckets = num_tuples_in_cache;

  fprintf(stderr, "Fanout: %u (%u bits)\n", fanout, FANOUT_BITS);

  //////////////////////////////////////
  // Allocate the hash table
  //////////////////////////////////////
  hashtable_t *hts[fanout];
  for (uint32_t i = 0; i < fanout; i++) {
    allocate_hashtable(&hts[i], nbuckets);
    hts[i]->skip_bits = FANOUT_BITS;
    hts[i]->hash_mask = (hts[i]->num_buckets - 1) << hts[i]->skip_bits;
    if (i == 0) {
      fprintf(stderr, "# buckets: %lu\n", hts[i]->num_buckets);
    }
  }

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

  arg_t arg = {.nthreads = nthreads, .tid = 0, .fanout = fanout, 
               .hts = hts, .build = relR, .probe = relS};

  //////////////////////////////////////
  // BUILD
  //////////////////////////////////////
  build_hashtable_st(&arg);

#ifdef DEBUG
  for (uint32_t i = 0; i < fanout; i++) {
    hashtable_t *ht = hts[i];
    uint32_t occ = 0;
    for (uint32_t j = 0; j < ht->num_buckets; j++) {
      occ += ht->buckets[j].flag;
    }
    fprintf(stderr, "HT%u load-factor: %.2lf \%\n", i, (double)occ/(double)ht->num_buckets);
  }
  //fprintf(stderr, "Num allocs: %llu\n", num_allocs);
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
  int64_t result = probe_hashtable_st(&arg);

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
               (relR->num_tuples+relS->num_tuples), result, &start, &end);
#endif


  // Cleanup
  for (uint32_t i = 0; i < fanout; i++) {
    hashtable_t *ht = hts[i];
    destroy_hashtable(ht);
  }
  //destroy_hashtable(ht);

  return result;
}
