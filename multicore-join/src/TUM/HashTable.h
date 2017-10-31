/*
 * Authors: Harald Lang (TUM), Viktor Leis (TUM)
 */
#ifndef _HASHTABLE_H_
#define _HASHTABLE_H_


//#define _MEMSET_HT_THREAD0_
//#define _MEMSET_HT_ALLTHREADS_
//#define _SPIN_LOCK_

#include <algorithm>
#include <assert.h>
#include <nmmintrin.h>
#include <pthread.h>
#include <sys/mman.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <iostream>

#include "Lock.h"
#include "Utils.h"
#include "types.h"
#include "Barrier.h"

using namespace std;

struct HTEntry {
    intkey_t hash;
    tuple_t t;
#if defined(_LOCK_IN_BUCKET_)
   Lock l;
#endif
};

#define COMMENT SLASH(/)
#define SLASH(s) /##s

#ifndef NO_TIMING
#define COUT cout
#else
#define COUT while(0) cout
#endif

class HashTable {
public:
    uint64_t mask;
    HTEntry *table;
    uint64_t count;
    bool lockFree;
    Lock *locks;
    size_t locksLength;
    size_t lockMask;
    size_t striping;
    uint64_t nBuckets;
	PThreadLockCVBarrier *barrier;

    HashTable(unsigned bits, bool lockFree, unsigned stripingBits, unsigned int nThreads) :
            table(NULL), count(0), lockFree(lockFree), locks(NULL) {
        nBuckets = 1ull << bits;
        striping = 1ull << stripingBits;
        locksLength = nBuckets / striping;
        lockMask = locksLength - 1;
        mask = nBuckets - 1;
		barrier = new PThreadLockCVBarrier(nThreads);

//		table = static_cast<HTEntry*>((void*) new char[nBuckets	* sizeof(HTEntry)]);
//		memset(table, 0, nBuckets * sizeof(HTEntry));
        COUT << "nBuckets=" << nBuckets << " mask=" << mask << endl;

#if defined(_IDHASH_)
        COUT << "ID hash" << endl;
#elif defined(_FIBHASH_)
        COUT << "FIB hash" << endl;
#elif defined(_CRCHASH_)
        COUT << "CRC hash" << endl;
#else
        COUT << "MURMUR hash" << endl;
#endif
    }

    ~HashTable() {
        munmap(table, (mask + 1) * sizeof(HTEntry));
        if (!lockFree) {
#if !defined(_LOCK_IN_BUCKET_)
            delete[] locks;
#endif
        }
    }

    void init(uint64_t threadId, uint64_t numThreads) {
#if defined(_MEMSET_HT_THREAD0_)
		if (threadId == 0) {
			COUT << "HT: initializing memory in thread 0" << endl;
		}

		if (threadId == 0) {
			table = static_cast<HTEntry*>(Utils::malloc_huge_and_set(nBuckets * sizeof(HTEntry)));
			if (!lockFree) {
				locks = new Lock[locksLength];
			}
		}
#elif defined(_MEMSET_HT_ALLTHREADS_)
		uint64_t htSize = nBuckets * sizeof(HTEntry);
		if (threadId == 0) {
			COUT << "HT: initializing memory in all threads" << endl;
			table = static_cast<HTEntry*>(Utils::malloc_huge(htSize));
		}
		barrier->Arrive();

		uint64_t memChunkSize = 1024*2;
		uint64_t numChunks = htSize / memChunkSize;

		intkey_t * initOrder = Utils::generateShuffledNumbers(numChunks, numChunks, 19650218ULL);
		for (uint64_t i = threadId; i < numChunks; i += numThreads) {
			memset(((char*)table) + (initOrder[i] * memChunkSize), 0, memChunkSize);
		}
		if (!lockFree) {
			locks = (Lock*)Utils::malloc_huge(locksLength * sizeof(Lock));
			numChunks = (locksLength * sizeof(Lock)) / memChunkSize;
			if (numChunks > 0) {
				initOrder = Utils::generateShuffledNumbers(numChunks, numChunks, 19790415ULL);
				for (uint64_t i = threadId; i < numChunks; i += numThreads) {
					memset(((char*)locks) + (initOrder[i] * memChunkSize), 0, memChunkSize);
				}
			}
			memset(locks, 0, locksLength * sizeof(Lock));
		}
		delete [] initOrder;
        barrier->Arrive();
#else
        if (threadId == 0) {
            COUT << "HT: dynamically initializing memory in all threads" << endl;
        }

        table = static_cast<HTEntry *>(Utils::malloc_huge(nBuckets * sizeof(HTEntry)));
        if (!lockFree) {
            locks = new Lock[locksLength];
        }
#endif

    }

#if defined(_IDHASH_)
	/** Identity Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return (k & mask);
	}
#elif defined(_FIBHASH_)
	/** Fibonacci Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return (k * 11400714819323198485ull) | (1ull<<((sizeof(intkey_t)*8-1)));
	}
#elif defined(_CRCHASH_)
	/** CRC Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return _mm_crc32_u64(0, k) | (1ull<<((sizeof(intkey_t)*8-1)));
	}
#else

    /** MurmurHash64A */
    inline intkey_t hashKey(intkey_t k) const {
        const intkey_t m = 0xc6a4a7935bd1e995;
        const int r = 47;
        intkey_t h = 0x8445d61a4e774912 ^(8 * m);
        k *= m;
        k ^= k >> r;
        k *= m;
        h ^= k;
        h *= m;
        h ^= h >> r;
        h *= m;
        h ^= h >> r;
        return h | (1ull << ((sizeof(intkey_t) * 8 - 1)));
    }

#endif


    inline Lock *getLock(uint64_t bucket) {
        Lock *lock;
#if defined(_LOCK_IN_BUCKET_)
		lock = &table[bucket / striping].l;
#else
        lock = &locks[bucket / striping];
#endif
        return lock;
    }

    inline bool tryWriteToBucket(uint64_t pos, intkey_t hash) {
        Lock *lock = getLock(pos);
        lock->lock();
        if (table[pos].hash == 0) {
            table[pos].hash = hash;
            lock->unlock();
            return true;
        } else {
            lock->unlock();
            return false;
        }
    }

    inline void insert(intkey_t key, value_t value) {
#if defined(_ATOMIC_)
		insertAtomic(key, value);
#elif defined(_NO_SYNC_)
		insertWithoutSync(key, value);
#else
        insertWithLock(key, value);
#endif
    }

    inline void insertWithLock(intkey_t key, value_t value) {
        intkey_t hash = hashKey(key);
        unsigned m = mask;
        uint64_t pos = hash & m;
        while (table[pos].hash || (!tryWriteToBucket(pos, hash))) {
            pos = (pos + 1) & m;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
    }

    inline void insertWithLockStriping(intkey_t key, value_t value) {
        intkey_t hash = hashKey(key);
        unsigned m = mask;
        uint64_t pos = hash & m;
        Lock *lock = getLock(pos);
        lock->lock();

        while (table[pos].hash) {
            pos = (pos + 1) & m;

            Lock *lock2 = getLock(pos);
            if (lock != lock2) {
                lock->unlock();
                lock = lock2;
                lock->lock();
            }

        }
        table[pos].hash = hash;
        table[pos].t.key = key;
        table[pos].t.payload = value;
        lock->unlock();
    }


    inline void insertAtomic(const intkey_t key, const value_t value) {
        const intkey_t hash = hashKey(key);
        const uint64_t m = mask;
        uint64_t pos = hash & m;
        while (table[pos].hash
                || (!__sync_bool_compare_and_swap(&table[pos].hash, 0, hash))) {
            pos = (pos + 1) & m;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
    }

    inline void insertWithoutSync(const intkey_t key, const value_t value) {
        intkey_t hash = hashKey(key);
        intkey_t m = mask;
        uint64_t pos = hash & m;
        while (table[pos].hash) {
            pos = (pos + 1) & m;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
        table[pos].hash = hash;
    }

    inline void insertPos(uint64_t pos, const intkey_t hash, const intkey_t key, const value_t value) {
#if defined(_ATOMIC_)
		insertPosAtomic(pos, hash, key, value);
#elif defined(_NO_SYNC_)
		insertPosWithoutSync(pos, hash, key, value);
#else
        insertPosWithLock(pos, hash, key, value);
#endif
    }

    inline void insertPosAtomic(uint64_t pos, const intkey_t hash, const intkey_t key, const value_t value) {
        while (table[pos].hash
                || (!__sync_bool_compare_and_swap(&table[pos].hash, 0, hash))) {
            pos = (pos + 1) & mask;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
    }

    inline void insertPosWithLock(uint64_t pos, const intkey_t hash, const intkey_t key, const value_t value) {
        while (table[pos].hash || (!tryWriteToBucket(pos, hash))) {
            pos = (pos + 1) & mask;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
    }

    inline void insertPosWithLockStriping(uint64_t pos, const intkey_t hash, const intkey_t key, const value_t value) {

        unsigned m = mask;

        Lock *lock = getLock(pos);
        lock->lock();

        while (table[pos].hash) {
            pos = (pos + 1) & m;

            Lock *lock2 = getLock(pos);
            if (lock != lock2) {
                lock->unlock();
                lock = lock2;
                lock->lock();
            }

        }

        table[pos].t.key = key;
        table[pos].t.payload = value;

        lock->unlock();
    }

    inline void insertPosWithoutSync(uint64_t pos, const intkey_t hash, const intkey_t key, const value_t value) {
        while (table[pos].hash) {
            pos = (pos + 1) & mask;
        }
        table[pos].t.key = key;
        table[pos].t.payload = value;
        table[pos].hash = hash;
    }

    inline void *lookup(const intkey_t key) const {
        uint64_t pos = hashKey(key) & mask;
        while ((table[pos].hash) && (table[pos].t.key != key)) {
            pos = (pos + 1) & mask;
        }
//		return table[pos].t.key == key ? &table[pos].t : NULL;
        return table[pos].hash ? &table[pos].t : NULL;
    }

    inline void *lookupPos(uint64_t pos, const intkey_t key) const {
        while ((table[pos].hash) && (table[pos].t.key != key)) {
            pos = (pos + 1) & mask;
        }
//		return table[pos].t.key == key ? &table[pos].t : NULL;
        return table[pos].hash ? &table[pos].t : NULL;
    }

    uint64_t chainLen(const intkey_t key) const {
        uint64_t i = 1;
        uint64_t pos = hashKey(key) & mask;
        while ((table[pos].hash) && (table[pos].t.key != key)) {
            pos = (pos + 1) & mask;
            i++;
        }
        return i;
    }


    inline void prefetchBuild(const intkey_t key) const {
        const uint64_t pos = hashKey(key) & mask;
        __builtin_prefetch(&table[pos], 1, 1);
    }

    inline void prefetchPosBuild(const uint64_t pos) const {
        __builtin_prefetch(&table[pos], 1, 1);
    }

    inline void prefetchProbe(const intkey_t key) const {
        const uint64_t pos = hashKey(key) & mask;
        __builtin_prefetch(&table[pos], 0, 1);
    }

    inline void prefetchPosProbe(const uint64_t pos) const {
        __builtin_prefetch(&table[pos], 0, 1);
    }

    inline void hash(const intkey_t key, intkey_t &hash, uint64_t &pos) const {
        hash = hashKey(key);
        pos = hash & mask;
    }

    inline void toPos(const intkey_t key, uint64_t &pos) const {
        pos = hashKey(key) & mask;
    }
};

#endif // _HASHTABLE_H_
