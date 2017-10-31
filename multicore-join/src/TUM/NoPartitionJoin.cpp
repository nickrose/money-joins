#include <algorithm>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/time.h>



//#define _IDHASH_ 1
#define _CRCHASH_ 1
#define _ATOMIC_ 1
#define _MEMSET_HT_ALLTHREADS_ 1

#include "Affinitizer.h"
#include "Barrier.h"
#include "HashTable.h"
#include "Lock.h"
#include "ProcessorMap.h"
#include "Utils.h"

#include "ParallelJoin.h"
#include "NoPartitionHashJoin.h"
#include "NoPartitionArrayJoin.h"
//#include "RadixHashJoin.h"
//#include "RadixHashJoinRandomScatter.h"
//#include "DMPSM.h"

using namespace std;

ProcessorMap processorMap;
Affinitizer affinitizer;

struct ThreadArgs {
	unsigned threadNum;
	ParallelJoin* join;
};

void* run(void *arg) {

	ThreadArgs* threadArgs = reinterpret_cast<ThreadArgs*>(arg);

	affinitizer.affinitize(threadArgs->threadNum);
	threadArgs->join->runInParallel(threadArgs->threadNum);

	delete threadArgs;
	return NULL;
}

template<bool is_checksum>
join_result_t
HYPER_NOPA(relation_t *R, relation_t *S, unsigned int threadCount)
{


	NoPartitionArrayJoin<is_checksum> * parallelJoin = new NoPartitionArrayJoin<is_checksum>(threadCount, R, S);

	// here we go
	pthread_t threads[1000];
	affinitizer.init(threadCount, processorMap.NumberOfProcessors());
	for (unsigned i = 0; i < threadCount; i++) {
		ThreadArgs* threadArgs = new ThreadArgs();
		threadArgs->join = parallelJoin;
		threadArgs->threadNum = i;
		pthread_create(&threads[i], NULL, run, reinterpret_cast<void*>(threadArgs));
	}
	for (unsigned i = 0; i < threadCount; i++) {
		pthread_join(threads[i], NULL);
	}
	join_result_t result = parallelJoin->getResult();
	delete parallelJoin;
	return result;
}

template join_result_t HYPER_NOPA<true>(relation_t *R, relation_t *S, unsigned int threadCount);
template join_result_t HYPER_NOPA<false>(relation_t *R, relation_t *S, unsigned int threadCount);

template<bool is_checksum>
/*join_result_t*/
int64_t
HYPER_NOP(relation_t *R, relation_t *S, int threadCount)
{


	NoPartitionHashJoin<is_checksum> * parallelJoin = new NoPartitionHashJoin<is_checksum>(threadCount, R, S);

	// here we go
	pthread_t threads[1000];
	affinitizer.init(threadCount, processorMap.NumberOfProcessors());
	for (unsigned i = 0; i < threadCount; i++) {
		ThreadArgs* threadArgs = new ThreadArgs();
		threadArgs->join = parallelJoin;
		threadArgs->threadNum = i;
		pthread_create(&threads[i], NULL, run, reinterpret_cast<void*>(threadArgs));
	}
	for (unsigned i = 0; i < threadCount; i++) {
		pthread_join(threads[i], NULL);
	}
	join_result_t result = parallelJoin->getResult();
	delete parallelJoin;
	//return result;
	return result.matches;
}

template /*join_result_t*/int64_t HYPER_NOP<true>(relation_t *R, relation_t *S, int threadCount);
template /*join_result_t*/int64_t HYPER_NOP<false>(relation_t *R, relation_t *S, int threadCount);
