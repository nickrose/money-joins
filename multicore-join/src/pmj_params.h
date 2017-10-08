#ifndef PMJ_PARAMS_H
#define PMJ_PARAMS_H
// DEBUG
// KEY_8B
// NO_TIMING
// PREFETCH_NPJ

/** Number of tuples that each bucket can hold */
#ifndef BUCKET_SIZE
#define BUCKET_SIZE 3
#endif

/** Size of system cache line in bytes */
#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

/** Pre-allocation size for overflow buffers */
#ifndef OVERFLOW_BUF_SIZE
#define OVERFLOW_BUF_SIZE 1024 
#endif

/** Should hashtable buckets be padded to cache line size */
#ifndef PADDED_BUCKET
#define PADDED_BUCKET 0 /* default case: not padded */
#endif

#ifndef PREFETCH_DISTANCE
#define PREFETCH_DISTANCE 16
#endif

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 1024
#endif

#endif /* PMJ_PARAMS_H */
