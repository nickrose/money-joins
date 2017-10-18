#include <inttypes.h>
#include <string>
#include <libcuckoo/cuckoohash_map.hh>

#include "common.h"

typedef cuckoohash_map<KeyType, ValueType, Hasher> hash_t;
typedef cuckoohash_map<std::string, ValueType, StringHasher> str_hash_t;

#define SETUP                                       \
  hash_t mt_hash;                                   \
  hash_t::locked_table hash = mt_hash.lock_table(); \
  str_hash_t mt_str_hash;                               \
  str_hash_t::locked_table str_hash = mt_str_hash.lock_table();

#define RESERVE_INT(size) hash.reserve(size)
#define RESERVE_STR(size) str_hash.reserve(size)
#define LOAD_FACTOR(map) map.load_factor()

#define INSERT_INT_INTO_HASH(key, value) \
  auto ret = hash.insert(key, value);    \
  if (!ret.second) {                     \
    printf("error on insert");           \
    exit(1);                             \
  }

#define DELETE_INT_FROM_HASH(key)
#define FIND_INT_EXISTING_FROM_HASH(key) \
  if (hash.find(key) == hash.end()) {    \
    printf("error");                     \
    exit(1);                             \
  }
#define FIND_INT_MISSING_FROM_HASH(key) \
  if (hash.find(key) != hash.end()) {   \
    printf("error");                    \
    exit(2);                            \
  }
#define FIND_INT_EXISTING_FROM_HASH_COUNT(key, count) \
  ValueType v;                                        \
  if (hash.find(key) != hash.end()) {                 \
    count++;                                          \
  }
#define CHECK_INT_ITERATOR_VALUE(iterator, value)

#define INSERT_STR_INTO_HASH(key, value) \
  auto ret = str_hash.insert(key, value);    \
  if (!ret.second) {                     \
    printf("error on insert");           \
    exit(1);                             \
  }

#define DELETE_STR_FROM_HASH(key)
#define FIND_STR_EXISTING_FROM_HASH(key) \
  if (str_hash.find(key) == str_hash.end()) {    \
    printf("error");                     \
    exit(1);                             \
  }
#define FIND_STR_MISSING_FROM_HASH(key) \
  if (str_hash.find(key) != str_hash.end()) {   \
    printf("error");                    \
    exit(2);                            \
  }
#define FIND_STR_EXISTING_FROM_HASH_COUNT(key, count) \
  ValueType v;                                        \
  if (str_hash.find(key) != str_hash.end()) {                 \
    count++;                                          \
  }

#include "template.c"
