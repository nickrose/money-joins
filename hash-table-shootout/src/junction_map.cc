#include <inttypes.h>
#include <string>
#include <junction/SingleMap_Leapfrog.h>

#include "common.h"

struct Int64KeyTraits {
  typedef size_t Hash;
  static const std::int64_t NullKey = -1;
  static const Hash NullHash = -1;    
  static inline Hash hash(std::int64_t key) {
  static Hasher hasher;
    return hasher(key);
  }
};

struct Int64ValueTraits {
  static const std::int64_t NullValue = -1;
};

typedef junction::SingleMap_Leapfrog<std::int64_t, const std::int64_t*, Int64KeyTraits, Int64ValueTraits> hash_t;
//typedef junction::SingleMap_Leapfrog<std::string, const std::int64_t*> str_hash_t;

#define SETUP hash_t hash; hash_t str_hash;

#define RESERVE_INT(size) 
#define RESERVE_STR(size) 
#define LOAD_FACTOR(map) map.load_factor()

#define INSERT_INT_INTO_HASH(key, value) hash.set(key, &value)
#define DELETE_INT_FROM_HASH(key) 
#define FIND_INT_EXISTING_FROM_HASH(key) if(hash.get(key) == nullptr) { printf("error"); exit(1); }
#define FIND_INT_MISSING_FROM_HASH(key) if(hash.get(key) != nullptr) { printf("error"); exit(2); }
#define FIND_INT_EXISTING_FROM_HASH_COUNT(key, count) \
            if(hash.get(key) != nullptr) { count++; }
#define CHECK_INT_ITERATOR_VALUE(iterator, value) 

#define INSERT_STR_INTO_HASH(key, value)
#define DELETE_STR_FROM_HASH(key) 
#define FIND_STR_EXISTING_FROM_HASH(key)
#define FIND_STR_MISSING_FROM_HASH(key) 
#define FIND_STR_EXISTING_FROM_HASH_COUNT(key, count) 

#include "template.c"
