#include <inttypes.h>
#include <string>
#include <unordered_map>
#include <tsl/array_map.h>

#include "common.h"

#include <experimental/string_view>
template<class CharT>
struct str_hash {
    std::size_t operator()(const CharT* key, std::size_t key_size) const {
        return std::hash<std::experimental::string_view>()(std::experimental::string_view(key, key_size));
    }
};

typedef std::unordered_map<int64_t, int64_t> hash_t;
typedef tsl::array_map<char, int64_t, str_hash<char>> str_hash_t;

#define SETUP hash_t hash; str_hash_t str_hash;

#define RESERVE_INT(size) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define RESERVE_STR(size) str_hash.reserve(size); 
#define LOAD_FACTOR(map) map.load_factor()

#define INSERT_INT_INTO_HASH(key, value) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define DELETE_INT_FROM_HASH(key) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define FIND_INT_EXISTING_FROM_HASH(key) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define FIND_INT_MISSING_FROM_HASH(key) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define FIND_INT_EXISTING_FROM_HASH_COUNT(key, count) { printf("array_map can't be used for integer benchmark"); exit(6); }
#define CHECK_INT_ITERATOR_VALUE(iterator, value) { printf("array_map can't be used for integer benchmark"); exit(6); }


#define INSERT_STR_INTO_HASH(key, value) str_hash.insert(key, value)
#define DELETE_STR_FROM_HASH(key) str_hash.erase(key);
#define FIND_STR_EXISTING_FROM_HASH(key) if(str_hash.find(key) == str_hash.end()) { printf("error"); exit(4); }
#define FIND_STR_MISSING_FROM_HASH(key) if(str_hash.find(key) != str_hash.end()) { printf("error"); exit(5); }
#define FIND_STR_EXISTING_FROM_HASH_COUNT(key, count) \
    if(str_hash.find(key) != str_hash.end()) { count++; }

#include "template.c"
