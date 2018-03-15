#ifndef HOPSCOTCH_JOIN_H
#define HOPSCOTCH_JOIN_H

#include "types.h" /* relation_t */

/** 
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 * 
 * @return number of result tuples
 */
int64_t HopscotchJoin(relation_t *relR, relation_t *relS, int nthreads);

#endif /* HOPSCOTCH_JOIN_H */
