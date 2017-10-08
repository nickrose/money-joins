#ifndef P_MONEY_JOIN_7_H
#define P_MONEY_JOIN_7_H

#include "types.h" /* relation_t */

/** 
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 * 
 * @return number of result tuples
 */
int64_t PMJ_7(relation_t *relR, relation_t *relS, int nthreads);

#endif /* P_MONEY_JOIN_7_H */
