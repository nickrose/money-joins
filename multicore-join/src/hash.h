#ifndef HASH_H
#define HASH_H

#ifdef MURMUR

static inline uint64_t Hash(uint64_t key) {
  //
  // The following code segment is copied from MurmurHash3, and is used
  // as an answer on the Internet:
  // http://stackoverflow.com/questions/5085915/what-is-the-best-hash-
  //   function-for-uint64-t-keys-ranging-from-0-to-its-max-value
  //
  key ^= key >> 33;
  key *= 0xff51afd7ed558ccd;
  key ^= key >> 33;
  key *= 0xc4ceb9fe1a85ec53;
  key ^= key >> 33;

  return key;
}
#else

//#warning "Using dumb hash"
// static inline uint64_t Hash(uint64_t key) {
// 	return key;
// }
#include <immintrin.h>

static inline uint64_t Hash(uint64_t key) {
  //return key;
  return _mm_crc32_u32(0x4c11db7, key);
  //return (key >> 7) ^ (key >> 13 ) ^ (key >> 21 ) ^ key;
}
#endif

#endif  // HASH_H
