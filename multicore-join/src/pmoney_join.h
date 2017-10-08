#ifndef P_MONEY_JOIN_H
#define P_MONEY_JOIN_H

#include "types.h" /* relation_t */

/** 
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 * 
 * @return number of result tuples
 */
int64_t PMJ(relation_t *relR, relation_t *relS, int nthreads);

#endif /* P_MONEY_JOIN_H */
