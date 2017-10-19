#include <inttypes.h>
#include <string>
#include <google/dense_hash_map>

#include "common.h"

typedef google::dense_hash_map<KeyType, ValueType, Hasher> hash_t;
typedef google::dense_hash_map<std::string, int64_t, StringHasher> str_hash_t;

#define SETUP hash_t hash; hash.max_load_factor(0.99f); \
			  hash.set_empty_key(-1); hash.set_deleted_key(-2); \
              str_hash_t str_hash; str_hash.max_load_factor(0.99f); \
              str_hash.set_empty_key(""); str_hash.set_deleted_key("d");

#define RESERVE_INT(size) hash.resize(size);
#define RESERVE_STR(size) str_hash.resize(size);
#define LOAD_FACTOR(map) map.load_factor()        

#define INSERT_INT_INTO_HASH(key, value) hash.insert(hash_t::value_type(key, value))
#define DELETE_INT_FROM_HASH(key) hash.erase(key)
#define FIND_INT_EXISTING_FROM_HASH(key) if(hash.find(key) == hash.end()) { printf("error"); exit(1); }
#define FIND_INT_MISSING_FROM_HASH(key) if(hash.find(key) != hash.end()) { printf("error"); exit(2); }
#define FIND_INT_EXISTING_FROM_HASH_COUNT(key, count) \
    if(hash.find(key) != hash.end()) { count++; }
#define CHECK_INT_ITERATOR_VALUE(iterator, value) if(iterator->second != value) { printf("error"); exit(3); }

#define INSERT_STR_INTO_HASH(key, value) str_hash.insert(str_hash_t::value_type(key, value))
#define DELETE_STR_FROM_HASH(key) str_hash.erase(key);
#define FIND_STR_EXISTING_FROM_HASH(key) if(str_hash.find(key) == str_hash.end()) { printf("error"); exit(4); }
#define FIND_STR_MISSING_FROM_HASH(key) if(str_hash.find(key) != str_hash.end()) { printf("error"); exit(5); }
#define FIND_STR_EXISTING_FROM_HASH_COUNT(key, count) \
    if(str_hash.find(key) != str_hash.end()) { count++; }
        
#include "template.c"
