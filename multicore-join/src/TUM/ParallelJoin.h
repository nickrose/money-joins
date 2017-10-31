/*
 * Authors: Harald Lang (TUM), Viktor Leis (TUM)
 */

#ifndef PARALLELJOIN_H_
#define PARALLELJOIN_H_

#include "Barrier.h"
#include "Lock.h"

class ParallelJoin {
protected:
	uint64_t threadCount;
	PThreadLockCVBarrier* barrier;
	Lock stopwatchMutex;

public:
	ParallelJoin(uint64_t threadCount) :
		threadCount(threadCount) {
		barrier = new PThreadLockCVBarrier(threadCount);
	}

	virtual ~ParallelJoin() {
		delete barrier;
	}

	virtual void runInParallel(uint64_t threadNum) = 0;

	void syncThreads() {
		barrier->Arrive();
	}
};

#endif /* PARALLELJOIN_H_ */
