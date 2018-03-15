
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

#include "hopscotch_join.h"
#include "hopscotch-map/src/hopscotch_map.h"

#include "pmj_params.h"         /* constant parameters */
#include "rdtsc.h"              /* startTimer, stopTimer */
#include "lock.h"               /* lock, unlock */
#include "cpu_mapping.h"        /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

#include "hash.h"
#include "generator.h"          /* numa_localize() */

#ifndef BARRIER_ARRIVE
/** barrier wait macro */
#define BARRIER_ARRIVE(B,RV)                            \
    RV = pthread_barrier_wait(B);                       \
    if(RV !=0 && RV != PTHREAD_BARRIER_SERIAL_THREAD){  \
        printf("Couldn't wait on barrier\n");           \
        exit(EXIT_FAILURE);                             \
    }
#endif

// An experimental feature to allocate input relations numa-local
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

struct Hasher {
  size_t operator()(const intkey_t &k) const { return Hash(k); }
};

using Map = tsl::hopscotch_map<intkey_t, value_t, Hasher>;

struct arg_t {
    int32_t             tid;
    Map           *     ht;
    relation_t          relR;
    relation_t          relS;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1, timer2, timer3;
    struct timeval start, end;
#endif
} ;

//===----------------------------------------------------------------------===//
// Allocate a hash table with the given number of buckets
//===----------------------------------------------------------------------===//
static inline std::unique_ptr<Map>
allocate_hashtable(uint32_t nelems) {
  auto *map = new Map();
  map->max_load_factor(0.5f);
  map->reserve(nelems);
  fprintf(stderr, "HT size: %.2lf MB\n",
          ((double)map->size_in_bytes() / 1024.0 / 1024.0));
  return std::unique_ptr<Map>{map};
}

//===----------------------------------------------------------------------===//
// Build a hash table using a single thread
//===----------------------------------------------------------------------===//
static inline void build_hashtable(Map *ht, relation_t *rel) {
  for(uint32_t i = 0; i < rel->num_tuples; i++) {
    ht->insert(Map::value_type(rel->tuples[i].key, rel->tuples[i].key));
  }
}

//===----------------------------------------------------------------------===//
// Probe the hash table using a single thread
//===----------------------------------------------------------------------===//
static inline int64_t probe_hashtable(Map *ht, relation_t *rel) {
  uint64_t matches = 0;
  for (uint32_t i = 0; i < rel->num_tuples; i++) {
    matches += (ht->find(rel->tuples[i].key) != ht->end() ? 1 : 0);
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
  fflush(stderr);
}

void *hopscotch_thread(void * param) {
  int rv;
  arg_t * args = (arg_t*) param;

#ifdef PERF_COUNTERS
  if(args->tid == 0){
      PCM_initPerformanceMonitor(NULL, NULL);
      PCM_start();
  }
#endif
    
  /* wait at a barrier until each thread starts and start timer */
  BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
  /* the first thread checkpoints the start time */
  if(args->tid == 0){
    gettimeofday(&args->start, NULL);
    startTimer(&args->timer1);
    startTimer(&args->timer2); 
    args->timer3 = 0; /* no partitionig phase */
  }
#endif

  /* insert tuples from the assigned part of relR to the ht */
  build_hashtable(args->ht, &args->relR);

  /* wait at a barrier until each thread completes build phase */
  BARRIER_ARRIVE(args->barrier, rv);
  
  if (args->tid == 0) {
    fprintf(stderr, "HT size after inserts: %.2lf MB\n",
            ((double)args->ht->size_in_bytes() / 1024.0 / 1024.0));
  }

#ifdef PERF_COUNTERS
  if(args->tid == 0){
    PCM_stop();
    PCM_log("========== Build phase profiling results ==========\n");
    PCM_printResults();
    PCM_start();
  }
  /* Just to make sure we get consistent performance numbers */
  BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
  /* build phase finished, thread-0 checkpoints the time */
  if(args->tid == 0){
    stopTimer(&args->timer2); 
  }
#endif

  /* probe for matching tuples from the assigned part of relS */
  args->num_results = probe_hashtable(args->ht, &args->relS);

#ifndef NO_TIMING
  /* for a reliable timing we have to wait until all finishes */
  BARRIER_ARRIVE(args->barrier, rv);

  /* probe phase finished, thread-0 checkpoints the time */
  if(args->tid == 0){
    stopTimer(&args->timer1); 
    gettimeofday(&args->end, NULL);
  }
#endif

#ifdef PERF_COUNTERS
  if(args->tid == 0) {
      PCM_stop();
      PCM_log("========== Probe phase profiling results ==========\n");
      PCM_printResults();
      PCM_log("===================================================\n");
      PCM_cleanup();
  }
  /* Just to make sure we get consistent performance numbers */
  BARRIER_ARRIVE(args->barrier, rv);
#endif

    return 0;
}

//===----------------------------------------------------------------------===//
// Run the algorithm
//===----------------------------------------------------------------------===//
int64_t HopscotchJoin(relation_t *relR, relation_t *relS, int nthreads) {
  int64_t result = 0;
  int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
  int i, rv;
  cpu_set_t set;
  arg_t args[nthreads];
  pthread_t tid[nthreads];
  pthread_attr_t attr;
  pthread_barrier_t barrier;
  
  uint32_t nbuckets = relR->num_tuples;
  auto ht = allocate_hashtable(nbuckets);
  
  numR = relR->num_tuples;
  numS = relS->num_tuples;
  numRthr = numR / nthreads;
  numSthr = numS / nthreads;

  rv = pthread_barrier_init(&barrier, NULL, nthreads);
  if(rv != 0){
      printf("Couldn't create the barrier\n");
      exit(EXIT_FAILURE);
  }

  pthread_attr_init(&attr);
  for(i = 0; i < nthreads; i++){
      int cpu_idx = get_cpu_id(i);

      CPU_ZERO(&set);
      CPU_SET(cpu_idx, &set);
      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

      args[i].tid = i;
      args[i].ht = ht.get();
      args[i].barrier = &barrier;

      /* assing part of the relR for next thread */
      args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
      args[i].relR.tuples = relR->tuples + numRthr * i;
      numR -= numRthr;

      /* assing part of the relS for next thread */
      args[i].relS.num_tuples = (i == (nthreads-1)) ? numS : numSthr;
      args[i].relS.tuples = relS->tuples + numSthr * i;
      numS -= numSthr;

      rv = pthread_create(&tid[i], &attr, hopscotch_thread, (void*)&args[i]);
      if (rv){
          printf("ERROR; return code from pthread_create() is %d\n", rv);
          exit(-1);
      }

  }

  for(i = 0; i < nthreads; i++){
      pthread_join(tid[i], NULL);
      /* sum up results */
      result += args[i].num_results;
  }

#ifndef NO_TIMING
  /* now print the timing results: */
  print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
               relR->num_tuples, relS->num_tuples, result,
               &args[0].start, &args[0].end);
#endif

  return result;
}
