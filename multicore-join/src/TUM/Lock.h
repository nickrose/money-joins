
#ifndef _SPINLOCK_H_
#define _SPINLOCK_H_

#if !defined(_SPIN_LOCK_) && !defined(_SPIN_LOCK_CAS_) && !defined(_PTHREAD_LOCK_)
#define _SPIN_LOCK_
#endif

#if defined(_SPIN_LOCK_)
class Lock {
	public:
		Lock() : mutex(0) { }

		/** Call blocks and returns only when it has the lock. */
		inline void lock()
		{
			while(tas(&mutex)) {
#if defined(__i386__) || defined(__x86_64__)
				__asm__ __volatile__ ("pause\n");
#endif
			}
		}

		/** Unlocks the lock object. */
		inline void unlock()
		{
			mutex = 0;
		}

	private:
		inline int tas(volatile char* lock)
		{
			register char res = 1;
#if defined(__i386__) || defined(__x86_64__)
			__asm__ __volatile__ (
					"lock xchgb %0, %1\n"
					: "+q"(res), "+m"(*lock)
					:
					: "memory", "cc");
#elif defined(__sparc__)
			__asm__ __volatile__ (
					"ldstub [%2], %0"
					: "=r"(res), "+m"(*lock)
					: "r"(lock)
					: "memory");
#else
#error TAS not defined for this architecture.
#endif
			return res;
		}

		volatile char mutex;
};
#endif

#if defined(_PTHREAD_LOCK_)
class Lock {
	public:
		Lock() {
			pthread_mutex_init(&mutex, NULL);
		}

		/** Call blocks and retunrs only when it has the lock. */
		inline void lock() {
			pthread_mutex_lock(&mutex);
		}

		/** Unlocks the lock object. */
		inline void unlock()
		{
			pthread_mutex_unlock(&mutex);
		}

	private:
		pthread_mutex_t mutex;
};
#endif

#if defined(_SPIN_LOCK_CAS_)
class Lock {
	public:
		Lock() : mutex(0) { }

		inline void lock() {
			while(! __sync_bool_compare_and_swap(&mutex, 0, 1)) {
#if defined(__i386__) || defined(__x86_64__)
				__asm__ __volatile__ ("pause\n");
#endif
			}
		}

		inline void unlock() {
			mutex = 0;
		}

	private:
		volatile char mutex;
};
#endif

#endif
