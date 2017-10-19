#include <inttypes.h>
#include <string>
#include <flat_hash_map.hpp>

#include "common.h"

struct power_of_two_hash : public Hasher {
    typedef ska::power_of_two_hash_policy hash_policy;
};

struct str_power_of_two_hash : public StringHasher {
    typedef ska::power_of_two_hash_policy hash_policy;
};

typedef ska::flat_hash_map<KeyType, ValueType, power_of_two_hash> hash_t;
typedef ska::flat_hash_map<std::string, int64_t, str_power_of_two_hash> str_hash_t;

#define SETUP hash_t hash; str_hash_t str_hash; hash.max_load_factor(0.99f); str_hash.max_load_factor(0.99f);

#define RESERVE_INT(size) hash.reserve(size);
#define RESERVE_STR(size) str_hash.reserve(size); 
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

