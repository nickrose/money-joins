#include <inttypes.h>
#include <string>
#include <libcuckoo/cuckoohash_map.hh>

#include "common.h"

typedef cuckoohash_map<KeyType, ValueType, Hasher> hash_t;
typedef cuckoohash_map<std::string, std::int64_t> str_hash_t;

#define SETUP  \
  hash_t hash; \
  hash_t str_hash;

#define RESERVE_INT(size) hash.reserve(size)
#define RESERVE_STR(size) str_hash.reserve(size)
#define LOAD_FACTOR(map) map.load_factor()

#define INSERT_INT_INTO_HASH(key, value) \
  if (!hash.insert(key, value)) {        \
    printf("error on insert");           \
    exit(1);                             \
  }
#define DELETE_INT_FROM_HASH(key)
#define FIND_INT_EXISTING_FROM_HASH(key) \
  ValueType v;                           \
  if (!hash.find(key, v)) {              \
    printf("error");                     \
    exit(1);                             \
  }
#define FIND_INT_MISSING_FROM_HASH(key) \
  ValueType v;                          \
  if (hash.find(key, v)) {              \
    printf("error");                    \
    exit(2);                            \
  }
#define FIND_INT_EXISTING_FROM_HASH_COUNT(key, count) \
  ValueType v;                                        \
  if (hash.find(key)) {                               \
    count++;                                          \
  }
#define CHECK_INT_ITERATOR_VALUE(iterator, value)

#define INSERT_STR_INTO_HASH(key, value)
#define DELETE_STR_FROM_HASH(key)
#define FIND_STR_EXISTING_FROM_HASH(key)
#define FIND_STR_MISSING_FROM_HASH(key)
#define FIND_STR_EXISTING_FROM_HASH_COUNT(key, count)

#include "template.c"
