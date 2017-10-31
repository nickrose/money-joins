#ifndef NOPARTITIONARRAYJOIN_H_
#define NOPARTITIONARRAYJOIN_H_

// increasing the key distance will cause a performance decrease with some hash functions (e.g. identity)
#ifndef KEY_DISTANCE
#define KEY_DISTANCE 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <types.h>

#include "HashTable.h"
#include "Utils.h"
#include "MersenneTwister.h"
#include "rdtsc.h"
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif
template<bool is_checksum>
class NoPartitionArrayJoin : public ParallelJoin {
private:
    tuple_t * R;
    uint64_t n;
    int ratioHoles;
    tuple_t * S;
    uint64_t m;

    uint64_t range;
    uint64_t k;
//    uint64_t maxKey;

    uint64_t* outputTupleCount;
    uint64_t* checkSums;
	uint64_t  time_usec;
	uint64_t  build_usec;
    int lockType;
    bool lockFree;
    int lockStripingBits;
    value_t* array;

//    uint64_t* buildKeys;
//    uint64_t* probeKeys;

    size_t prefetchOffset;
    struct PrefetchCache {
        uint64_t hash;
        uint64_t pos;
    };

public:
    NoPartitionArrayJoin(uint64_t threadCount, relation_t * relR, relation_t * relS) :
            ParallelJoin(threadCount),
            R(relR->tuples), n(relR->num_tuples), ratioHoles(/*relR->ratio_holes*/0), S(relS->tuples), m(relS->num_tuples), prefetchOffset(16) {

		k = 0;
        lockType = 0;
        lockFree = true;
        lockStripingBits = 0;


#if defined(_ATOMIC_)
		COUT << "Using lock free implementation (atomic CAS)..." << endl;
#elif defined(_NO_SYNC_)
		COUT << "Using implementation without any synchronization... (THIS PRODUCES WRONG RESULTS! Used just for comparison)" << endl;
		lockType = -1;
#else
        lockFree = false;

#if defined(_SPIN_LOCK_)
        COUT << "Using spin-locks (TAS)..." << endl;
		lockType = 1;
#elif defined(_SPIN_LOCK_CAS_)
		COUT << "Using spin-locks (CAS)..." << endl;
		lockType = 2;
#elif defined(_PTHREAD_LOCK_)
		COUT << "Using pthread-locks..." << endl;
		lockType = 3;
#endif


#if defined(_LOCK_IN_BUCKET_)
		COUT << "Storing locks directly within the buckets..." << endl;
		lockType += 10;
#else
        COUT << "Storing locks in a separate array..." << endl;
#endif

        COUT << "lock size = " << sizeof(Lock) << " byte(s)" << endl;

      //  if (argCount > 3) {
      //      lockStripingBits = atoi(args[3]);
      //  }
#endif

#if defined(_PREFETCH_)
	/*	if (lockFree) {
			if (argCount > 3) {
				prefetchOffset = atoi(args[3]);
			}
		} else {
			if (argCount > 4) {
				prefetchOffset = atoi(args[4]);
			}
		}*/
		prefetchOffset = Utils::nextPowerOfTwo(prefetchOffset);
		cout << "Using prefetching with offset = " << prefetchOffset << endl;
#endif


        outputTupleCount = new uint64_t[threadCount];
		if (is_checksum)
			checkSums = new uint64_t[threadCount];
        for (uint64_t i = 0; i < threadCount; i++) {
            outputTupleCount[i] = 0;
			if (is_checksum)
				checkSums[i] = 0;
        }
        range = n * ratioHoles;

    }

    ~NoPartitionArrayJoin() {
        delete [] outputTupleCount;
		delete [] checkSums;
        delete [] array;
    }

    join_result_t getResult() {
        join_result_t ret = {0,0, time_usec, build_usec, time_usec - build_usec};
        for (unsigned int i = 0; i < threadCount; ++i) {
            ret.matches += outputTupleCount[i];
			if (is_checksum)
				ret.checksum += checkSums[i];
        }
        return ret;
    }

    void runInParallel(uint64_t threadNum) {
        // initialize hash table
        {
			if (threadNum == 0)
			{
				array = static_cast<value_t*>(malloc(sizeof(value_t) * (range+1)));
				if (array == NULL)
				{
					printf("[ERROR] array allocation failed\n");
					exit(1);
				}
			}
            syncThreads();
			uint64_t numPerThread = (range+1) / threadCount;
			uint64_t myPartNum = (threadNum == threadCount - 1) ? (range+1 - threadNum*numPerThread) : numPerThread;
			uint64_t myOffset = threadNum * numPerThread;
			memset(array + myOffset, 0, sizeof(value_t) * myPartNum);
        }

        syncThreads();

        struct timeval start_timeval, end_timeval, build_timeval;
#ifdef ALGO_TIME
        if (threadNum == 0)
			gettimeofday(&start_timeval, NULL);
#endif
#ifdef PERF_COUNTERS
        if (threadNum == 0){
                PCM_initPerformanceMonitor(NULL, NULL);
	            PCM_start();
            }
#endif
        // build phase ------------------------------------------------------------------

#ifndef NO_TIMING
        double durationBuildPhase = 0.0;
        uint64_t cyclesBuildPhase = 0;
        double buildPerformance = 0.0;
#endif
        {

#ifndef NO_TIMING
            double start;
            if (threadNum == 0) {
                start = Utils::gettime();
                startTimer(&cyclesBuildPhase);
            }
#endif

            const uint64_t per_thread = n / threadCount;
            const uint64_t startIndex = threadNum * per_thread;
            const uint64_t endIndex = (threadNum == threadCount - 1) ? n : (threadNum + 1) * per_thread;

/*
#if defined(_PREFETCH_CACHE_)
			PrefetchCache* cache = new PrefetchCache[prefetchOffset];
			const uint64_t cacheMask = prefetchOffset - 1;

			// fill cache
			for (uint32_t i = startIndex; i < startIndex + prefetchOffset; i++) {
				const uint64_t cacheIndex = i & cacheMask;
				ht->hash(R[i].key, cache[cacheIndex].hash, cache[cacheIndex].pos);
				ht->prefetchPosBuild(cache[cacheIndex].pos);
			}

			// build - [startIndex ... endIndex - prefetchOffset)
			uint64_t prefetchIndex = startIndex + prefetchOffset; 	// adopted from the ETH NPO join algorithm
			const uint64_t endIndexPrefetch = endIndex - prefetchOffset;
			for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {
				// re-use hash values and positions, which have been already computed
				const uint64_t cacheIndex = i & cacheMask;
				ht->insertPos(cache[cacheIndex].pos, cache[cacheIndex].hash, R[i].key, R[i].payload);

				// prefetching
				ht->hash(R[prefetchIndex].key, cache[cacheIndex].hash, cache[cacheIndex].pos);
				ht->prefetchPosBuild(cache[cacheIndex].pos);
				prefetchIndex++;
			}
			for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {
				const uint64_t cacheIndex = i & cacheMask;
				ht->insertPos(cache[cacheIndex].pos, cache[cacheIndex].hash, R[i].key, R[i].payload);
			}

			delete [] cache;
#elif defined(_PREFETCH_)
			uint64_t prefetchIndex = startIndex + prefetchOffset; 	// adopted from the ETH NPO join algorithm
			for (uint64_t i = startIndex; i < endIndex; i++) {
				if (prefetchIndex < endIndex) {
					ht->prefetchBuild(R[prefetchIndex++].key);
				}
				ht->insert(R[i].key, R[i].payload);
			}
#else*/
            for (uint64_t i = startIndex; i < endIndex; i++) {
                array[R[i].key] = R[i].payload;
            }
//#endif

            syncThreads();
#ifdef ALGO_TIME
        if (threadNum == 0)
			gettimeofday(&build_timeval, NULL);
#endif

#ifndef NO_TIMING
            if (threadNum == 0) {
                cout << "build finished" << endl;
                stopTimer(&cyclesBuildPhase);
                durationBuildPhase = Utils::gettime() - start;
            }
#endif
#ifdef PERF_COUNTERS
            if(threadNum == 0){
        PCM_stop();
        PCM_log("======= Partitioning phase profiling results ======\n");
        PCM_printResults();
        PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
#endif
        }
        // full barrier
        __sync_synchronize();

        // probe phase ------------------------------------------------------------------
#ifndef NO_TIMING
        double durationProbePhase = 0.0;
        uint64_t cyclesProbePhase = 0;
        double probePerformance;
#endif
        {

#ifndef NO_TIMING
            double start;
            if (threadNum == 0) {
                start = Utils::gettime();
                startTimer(&cyclesProbePhase);
            }
#endif
            const uint64_t per_thread = m / threadCount;
            const uint64_t startIndex = threadNum * per_thread;
            const uint64_t endIndex = threadNum == threadCount - 1 ? m : (threadNum + 1) * per_thread;
            uint64_t count = 0;
            uint64_t checksum = 0;

			/*
#if defined(_PREFETCH_CACHE_)
			PrefetchCache* cache = new PrefetchCache[prefetchOffset];
			const uint64_t cacheMask = (1 << Utils::unsignedLog2(prefetchOffset)) - 1;

			// fill cache
			for (uint32_t i = startIndex; i < startIndex + prefetchOffset; i++) {
				const uint64_t cacheIndex = i & cacheMask;
				ht->hash(S[i].key, cache[cacheIndex].hash, cache[cacheIndex].pos);
				ht->prefetchPosProbe(cache[cacheIndex].pos);
			}

			// probe [0 .. n - prefetchOffset)
			uint64_t prefetchIndex = startIndex + prefetchOffset; 	// adopted from the ETH NPO join algorithm
			const uint64_t endIndexPrefetch = endIndex - prefetchOffset;
			for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {
				// re-use hash values and position, which have been already computed during prefetching
				const uint64_t cacheIndex = i & cacheMask;
				tuple_t *found_tuple = NULL;
				if ((found_tuple = reinterpret_cast<tuple_t*>(ht->lookupPos(cache[cacheIndex].pos, S[i].key))) != NULL) {
					count++;
					if (is_checksum)
						checksum += found_tuple->payload + S[i].payload;
				}

				// prefetching
				ht->hash(S[prefetchIndex].key, cache[cacheIndex].hash, cache[cacheIndex].pos);
				ht->prefetchPosProbe(cache[cacheIndex].pos);
				prefetchIndex++;
			}
			for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {
				const uint64_t cacheIndex = i & cacheMask;
				tuple_t *found_tuple = NULL;
				if ((found_tuple = reinterpret_cast<tuple_t*>(ht->lookupPos(cache[cacheIndex].pos, S[i].key))) != NULL) {
					count++;
					if (is_checksum)
						checksum += found_tuple->payload + S[i].payload;
				}
			}

			delete [] cache;

#elif defined(_PREFETCH_)
			size_t prefetchIndex = startIndex + prefetchOffset; // adopted from the ETH NPO join algorithm
			for (uint64_t i = startIndex; i < endIndex;	i++) {
				if (prefetchIndex < endIndex) {
					ht->prefetchProbe(S[prefetchIndex++].key);
				}
				tuple_t *found_tuple = NULL;
				if ((found_tuple = reinterpret_cast<tuple_t*>(ht->lookup(S[i].key))) != NULL) {
					count++;
					if (is_checksum)
						checksum += found_tuple->payload + S[i].payload;
				}
			}
#else*/
            for (uint64_t i = startIndex; i < endIndex;	i++) {
				value_t found_value = 0;
                if (S[i].key > 0 && S[i].key <= range && (found_value = array[S[i].key]) != 0) {
#if defined(DEBUG)
					tuple_t * t = (tuple_t *)ht->lookup(S[i].key);
					stringstream s;
					s << "(" << S[i].key << ", " << t->payload << ", " << S[i].key << ")" << endl;
					cout << s.str();
#endif

                    count++;
					if (is_checksum)
						checksum += found_value + S[i].payload;
                }
            }
//#endif
            outputTupleCount[threadNum] = count;
			if (is_checksum)
				checkSums[threadNum] = checksum;

#ifdef ALGO_TIME
            syncThreads();
			if (threadNum == 0) {
				gettimeofday(&end_timeval, NULL);
				time_usec = diff_usec(&start_timeval, &end_timeval);
				build_usec = diff_usec(&start_timeval, &build_timeval);
			}
#endif
#ifndef NO_TIMING
            syncThreads();
            if (threadNum == 0) {
                stopTimer(&cyclesProbePhase);
                durationProbePhase = Utils::gettime() - start;
            }
#endif
#ifdef PERF_COUNTERS
            if(threadNum == 0) { PCM_stop();
        PCM_log("=========== Build+Probe profiling results =========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
#endif
        }

#ifndef NO_TIMING
        // print report ------------------------------------------------------------------
        {
            if (threadNum == 0) {
                double durationTotal = durationBuildPhase + durationProbePhase;

                buildPerformance = ((n / (1024.0 * 1024.0)) / durationBuildPhase);
                probePerformance = ((m / (1024.0 * 1024.0)) / durationProbePhase);
                double overallPerformance = (((m + n) / (1024.0 * 1024.0)) / durationTotal);
                cout << "overall  => " << overallPerformance
                        << " M tuples per second." << endl;

                cout << "n=" << n << " m=" << m << " k=" << k << " threads=" << threadCount
                        << " joinKeyDistance=" << KEY_DISTANCE
                        << " build=" << buildPerformance
                        << " probe=" << probePerformance
                        << " total=" << overallPerformance
                        << endl;

                cout << "Performance (M tuples per second): ";

                cout << "part: " << 0;
                cout << " build: " << buildPerformance;
                cout << " probe: " << probePerformance;
                cout << " overall:" << overallPerformance;
                cout << endl << endl;

                uint64_t numOutputTuples = 0;
                for (unsigned int i = 0; i < threadCount; i++) {
                    numOutputTuples += outputTupleCount[i];
                }

                cout << "number of output tuples: " << numOutputTuples << endl;
                cout << "cycles per output tuple: " << ((cyclesBuildPhase + cyclesProbePhase) / numOutputTuples) << endl;

                cout << endl << "RUNTIME TOTAL, BUILD+PART, PART (cycles): " << endl;
                cout << "=cycles= " << (cyclesBuildPhase + cyclesProbePhase)
                        << "\t" << cyclesBuildPhase << "\t" << cyclesProbePhase << endl;

                cout << "Timers: TOTAL\tPART\tBUILD\tPROBE (seconds)" << endl;
                cout << "=wallclock= " << durationTotal << "\t"
                        << 0 << "\t"
                        << durationBuildPhase << "\t"
                        << durationProbePhase << endl;

                cout << "Done." << endl;
            }
        }
#endif
    }

};


#endif
