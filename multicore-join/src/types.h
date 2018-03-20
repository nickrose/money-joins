/**
 * @file    types.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Tue May 22 16:43:30 2012
 * @version $Id: types.h 3017 2012-12-07 10:56:20Z bcagri $
 * 
 * @brief  Provides general type definitions used by all join algorithms.
 * 
 * 
 */
#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <stdio.h>              /* printf */

/**
 * @defgroup Types Common Types
 * Common type definitions used by all join implementations.
 * @{
 */

#ifdef KEY_8B /* 64-bit key/value, 16B tuples */
typedef int64_t intkey_t;
typedef int64_t value_t;
#else /* 32-bit key/value, 8B tuples */
typedef int32_t intkey_t;
typedef int32_t value_t;
#endif

typedef struct tuple_t    tuple_t;
typedef struct relation_t relation_t;

/** Type definition for a tuple, depending on KEY_8B a tuple can be 16B or 8B */
struct tuple_t {
    intkey_t key;
    value_t  payload;
};

/**
 * Type definition for a relation. 
 * It consists of an array of tuples and a size of the relation.
 */
struct relation_t {
  tuple_t * tuples;
  uint32_t  num_tuples;
};

struct join_result_t {
  uint64_t matches;
  uint64_t checksum;
  uint64_t time_usec;
  uint64_t part_usec;
  uint64_t join_usec;
};

/** @} */

/** print out the execution time statistics of the join */
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
  double build_cpt = (double)build_cycles / (double)num_build;
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

#endif /* TYPES_H */
