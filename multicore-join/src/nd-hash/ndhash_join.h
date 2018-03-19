#ifndef NDHASH_JOIN_H
#define NDHASH_JOIN_H

#include "types.h" /* relation_t */

/** 
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 * 
 * @return number of result tuples
 */
int64_t NdHashJoin(relation_t *relR, relation_t *relS, int nthreads);

#endif /* NDHASH_JOIN_H */
