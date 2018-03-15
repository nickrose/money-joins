#ifndef HASH_H
#define HASH_H

#undef MURMUR
#undef CRCHASH
#undef FIBHASH
#undef IDHASH

#define MULTHASH

#if defined(MURMUR)
// MURMUR
#warning "Using Murmur"
static inline uint64_t Hash(intkey_t k) {
  #if 0
  k ^= k >> 16;
  k *= 0x85ebca6b;
  k ^= k >> 13;
  k *= 0xc2b2ae35;
  k ^= k >> 16;
  return k;
  #endif
  #if 1
  uint64_t kk = (uint64_t)k;
  //
  // The following code segment is copied from MurmurHash3, and is used
  // as an answer on the Internet:
  // http://stackoverflow.com/questions/5085915/what-is-the-best-hash-
  //   function-for-uint64-t-keys-ranging-from-0-to-its-max-value
  //
  kk ^= kk >> 33;
  kk *= 0xff51afd7ed558ccd;
  kk ^= kk >> 33;
  kk *= 0xc4ceb9fe1a85ec53;
  kk ^= kk >> 33;
  return kk;
  #endif

  #if 0
  const uint64_t m = 0xc6a4a7935bd1e995;
  const int r = 47;
  uint64_t h = 0x8445d61a4e774912 ^(8 * m);
  k *= m;
  k ^= k >> r;
  k *= m;
  h ^= k;
  h *= m;
  h ^= h >> r;
  h *= m;
  h ^= h >> r;
  return h | (1ull << ((sizeof(uint64_t) * 8 - 1)));
  #endif
}
#elif defined(CRCHASH)
// CRC
#warning "Using CRC"
#include <immintrin.h>
static inline intkey_t Hash(const intkey_t k) {
  return _mm_crc32_u64(0, k) | (1ull<<((sizeof(intkey_t)*8-1)));
}
#elif defined(FIBHASH)
// FIB
#warning "Using Fib Hash"
static inline intkey_t Hash(intkey_t k) {
  return (k * 11400714819323198485ull) | (1ull<<((sizeof(intkey_t)*8-1)));
}
#elif defined(MULTHASH)
// Multiplictative
static inline intkey_t Hash(const intkey_t k) {
  //return k * 2654435761LU; //golden ratio for 32 bits constant
  return k * 11400714819323198485ULL; //golden ratio for 64 bits constant
}
#elif defined(IDHASH)
/// ID HASH
#warning "Using ID Hash"
static inline uint64_t Hash(uint64_t key) {
  return key;
}
#else
/// XOR
#warning "Using XOR+Shift Hash"
static inline uint64_t Hash(uint64_t key) {
  return (key >> 7) ^ (key >> 13 ) ^ (key >> 21 ) ^ key;
}
#endif

#endif  // HASH_H
