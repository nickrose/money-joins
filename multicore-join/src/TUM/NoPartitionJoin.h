#ifndef HYPER_NOPARTITIONJOIN_H_
#define HYPER_NOPARTITIONJOIN_H_

#include "types.h"

template<bool is_checksum>
/*join_result_t*/int64_t HYPER_NOP(relation_t *, relation_t *, int);


template<bool is_checksum>
join_result_t HYPER_NOPA(relation_t *, relation_t *, unsigned int);



#endif
