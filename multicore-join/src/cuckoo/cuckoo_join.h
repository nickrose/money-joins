#ifndef CUCKOO_JOIN_H
#define CUCKOO_JOIN_H

#include "types.h" /* relation_t */

/** 
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 * 
 * @return number of result tuples
 */
int64_t CuckooJoin(relation_t *relR, relation_t *relS, int nthreads);

#endif /* CUCKOO_JOIN_H */
