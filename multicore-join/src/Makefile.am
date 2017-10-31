AUTOMAKE_OPTIONS = foreign

CC=g++

LANG =
DEFINES = 
MYLDFLAGS =
MISC =  -I/home/pmenon/work/multicore-joins/sniper-6.1/include

if KEY_8B
DEFINES += -DKEY_8B
endif

if PERF_COUNTERS
DEFINES += -DPERF_COUNTERS
MYLDFLAGS += -lperf
#MISC += -I../lib/intel-pcm-1.7/ -L../lib/intel-pcm-1.7/lib/
MISC += -I../lib/pcm/ -L../lib/pcm/lib
LIBS += -lperf
LANG += -std=c++11
else
#LANG += -std=c99
LANG += -std=c++11
endif

if DEBUG
DEFINES += -DDEBUG
endif

if PADDEDBUCKET
DEFINES += -DPADDED_BUCKET=1
else
DEFINES += -DPADDED_BUCKET=0
endif

if NOTIMING
DEFINES += -DNO_TIMING
endif

if SYNCSTATS
DEFINES += -DSYNCSTATS
endif

if SKEW_HANDLING
DEFINES += -DSKEW_HANDLING
endif

if PREFETCH_NPJ
DEFINES += -DPREFETCH_NPJ
endif

if USE_SWWC_OPTIMIZED_PART
DEFINES += -DUSE_SWWC_OPTIMIZED_PART
endif

if SEPARATE_CHAIN
DEFINES += -DSEPARATE_CHAIN
endif

# if available, add the following flags
# -mavx, -mtune=niagara2, -mcpu=ultrasparc
AM_CFLAGS = -Wall -ansi -Wno-format -pedantic \
	    -g -msse4.2 -mavx -mavx2 -march=native $(LANG) $(DEFINES)
AM_CXXFLAGS = $(AM_CFLAGS)

AM_LDFLAGS =
CFLAGS = -O3 -g -ggdb
CXXFLAGS = -O3
bin_PROGRAMS = mchashjoins
mchashjoins_SOURCES =	pmoney_join_9.h pmoney_join_9.c \
			TUM/Affinitizer.cpp TUM/Barrier.cpp TUM/NoPartitionJoin.cpp TUM/ProcessorMap.cpp \
			pmoney_join_8.h pmoney_join_8.c \
			hyper_join.h hyper_join.c \
			pmoney_join_7.h pmoney_join_7.c \
			pmoney_join_6.h pmoney_join_6.c \
			pmoney_join_5.h pmoney_join_5.c \
                        pmoney_join_4.h pmoney_join_4.c \
                        pmoney_join_3.h pmoney_join_3.c \
			pmoney_join_2.h pmoney_join_2.c \
			pmj_params.h pmoney_join.h pmoney_join.c \
 			npj_params.h prj_params.h types.h npj_types.h 	\
			no_partitioning_join.h no_partitioning_join.c 	\
			parallel_radix_join.h parallel_radix_join.c   	\
			perf_counters.h perf_counters.c			\
			cpu_mapping.h cpu_mapping.c 			\
			genzipf.h genzipf.c generator.h generator.c 	\
			lock.h rdtsc.h task_queue.h barrier.h affinity.h\
			main.c

mchashjoins_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS) -pthread
mchashjoins_CFLAGS  = $(AM_CFLAGS) $(MISC)
