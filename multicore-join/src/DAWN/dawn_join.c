
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>              /* CPU_ZERO, CPU_SET */
#include <pthread.h>            /* pthread_* */
#include <string.h>             /* memset */
#include <stdio.h>              /* printf */
#include <stdlib.h>             /* memalign */
#include <sys/time.h>           /* gettimeofday */
#include <memory>

#include "dawn_join.h"
#include "CuckooMap.h"
#include "pmj_params.h"         /* constant parameters */
#include "rdtsc.h"              /* startTimer, stopTimer */
#include "lock.h"               /* lock, unlock */
#include "cpu_mapping.h"        /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

#include "hash.h"
#include "generator.h"          /* numa_localize() */

#define PREFETCH

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline std::unique_ptr<CuckooMap<value_t>> 
allocate_hashtable(uint32_t nelems) {
  uint32_t cap = (((double)nelems) / 0.85);
  return std::unique_ptr<CuckooMap<value_t>>{new CuckooMap<value_t>(cap)};
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//
static inline void build_hashtable_st(CuckooMap<value_t> &ht, relation_t *rel) {
  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    ht.insert(rel->tuples[i].key, rel->tuples[i].key);
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable_st(CuckooMap<value_t> &ht, relation_t *rel) {
  uint64_t matches = 0;
#ifdef PREFETCH
  uint32_t pf_idx = 64;
#endif
  for (uint32_t i = 0; i < rel->num_tuples; i++) {
#ifdef PREFETCH
    if (pf_idx < rel->num_tuples) {
      ht.prefetch_buckets(rel->tuples[pf_idx++].key);
    }
#endif
    matches += (ht.get(rel->tuples[i].key).value == rel->tuples[i].key ? 1 : 0);
  }
  return matches;
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t DawnJoin(relation_t *relR, relation_t *relS, int nthreads) {

#ifndef NO_TIMING
  struct timeval start, end;
  uint64_t partition_time, build_time, probe_time;
#endif

  uint32_t nbuckets = relR->num_tuples;
  auto ht = allocate_hashtable(nbuckets);

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
  build_hashtable_st(*ht, relR);

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
  int64_t result = probe_hashtable_st(*ht, relS);

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

  return result;
}
