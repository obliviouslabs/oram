#pragma once
#include <immintrin.h>
#include <inttypes.h>

#include <cstring>
#include <typeinfo>
#include <utility>

#include "cpp_extended.hpp"
INLINE void CSWAP8(const uint64_t cond, uint64_t& guy1, uint64_t& guy2) {
  asm volatile(
      "test %[mcond], %[mcond]\n\t"
      "mov %[i1], %%r9\n\t"
      "cmovnz %[i2], %[i1]\n\t"
      "cmovnz %%r9, %[i2]\n\t"
      : [i1] "=r"(guy1), [i2] "=r"(guy2)
      : [mcond] "r"(cond), "[i1]"(guy1), "[i2]"(guy2)
      : "r9");
}

INLINE void CMOV8_internal(const uint64_t cond, uint64_t& guy1,
                           const uint64_t& guy2) {
  asm volatile(
      "test %[mcond], %[mcond]\n\t"
      "cmovnz %[i2], %[i1]\n\t"
      : [i1] "=r"(guy1)
      : [mcond] "r"(cond), "[i1]"(guy1), [i2] "r"(guy2)
      :);
}

INLINE void CMOV4_internal(const uint64_t cond, uint32_t& guy1,
                           const uint32_t& guy2) {
  asm volatile(
      "test %[mcond], %[mcond]\n\t"
      "cmovnz %[i2], %[i1]\n\t"
      : [i1] "=r"(guy1)
      : [mcond] "r"(cond), "[i1]"(guy1), [i2] "r"(guy2)
      :);
}

INLINE void CMOV1(const bool& cond, uint8_t& val1, const uint8_t& val2) {
  uint32_t r1 = 0 | val1;
  uint32_t r2 = 0 | val2;
  CMOV4_internal(cond, r1, r2);
  val1 = r1 & 0xff;
}

INLINE void CMOV2(const bool& cond, uint16_t& val1, const uint16_t& val2) {
  uint32_t r1 = 0 | val1;
  uint32_t r2 = 0 | val2;
  CMOV4_internal(cond, r1, r2);
  val1 = r1 & 0xffff;
}

INLINE void CMOV4(const bool& cond, uint32_t& val1, const uint32_t& val2) {
  CMOV4_internal(cond, val1, val2);
}

INLINE void CMOV8(const bool& cond, uint64_t& val1, const uint64_t& val2) {
  CMOV8_internal(cond, val1, val2);
}

INLINE void CMOV_BOOL(const uint64_t& cond, bool& val1, const bool& val2) {
  uint32_t v1 = val1;
  CMOV4_internal(cond, v1, val2);
  val1 = v1;
}

// Note: Possible generic implementation of CMOV, that we don't use.
// template<typename T>
// void ObliMov(bool mov, T* guy1, T* guy2) {
//   static_assert(sizeof(T)%8 == 0);
//   uint64_t* curr1 = (uint64_t*)guy1;
//   uint64_t* curr2 = (uint64_t*)guy2;
//   for (uint64_t i = 0; i < sizeof(T) / 8; ++i) {
//     CMOV(mov, *curr1, *curr2);
//     curr1++;
//     curr2++;
//   }
// }

template <typename T>
INLINE void CMOV(const uint64_t& cond, T& val1, const T& val2) {
  Assert(false,
         "This should mov not be compiled check that you implemented CMOV"
         "and that if you used overloading that you called "
         "OVERLOAD_TSET_CXCHG. For type: ",
         typeid(T).name());
  if (cond) {
    val1 = val2;
  }
}

template <typename T>
INLINE void TSET(bool selector, T& A, const T& B, const T& C) {
  CMOV(selector, A, C);
  CMOV(!selector, A, B);
}

template <typename T>
INLINE void CTSET(bool condition, bool selector, T& A, const T& B, const T& C) {
  CMOV(condition * selector, A, C);
  CMOV(condition * !selector, A, B);
}

template <typename T>
INLINE void CXCHG(const uint64_t& cond, T& A, T& B) {
  const T C = A;
  CMOV(cond, A, B);
  CMOV(cond, B, C);
}

template <const uint64_t sz>
INLINE void CXCHG_internal(const bool cond, void* vec1, void* vec2) {
  static_assert(sz <= 64);
#if defined(__AVX512VL__)
  const __mmask8 blend_mask = (__mmask8)(!cond) - 1;
#endif
  if constexpr (sz == 64) {
#if false && defined(__AVX512VL__)
    /* alternative implementation
    __m512i vec1_temp, vec2_temp;
    std::memcpy(&vec1_temp, vec1, 64);
    std::memcpy(&vec2_temp, vec2, 64);
    const __m512i& vec1_after_swap =
        _mm512_mask_blend_epi64(blend_mask, vec1_temp, vec2_temp);
    const __m512i& vec2_after_swap =
        _mm512_mask_blend_epi64(blend_mask, vec2_temp, vec1_temp);
    std::memcpy(vec1, &vec1_after_swap, 64);
    std::memcpy(vec2, &vec2_after_swap, 64);
    */
    __m512i vec1_temp, vec2_temp;
    __m512i temp;
    std::memcpy(&vec1_temp, vec1, 64);
    std::memcpy(&vec2_temp, vec2, 64);
    __m512i mask = _mm512_set1_epi32(-cond);  // Create a mask based on the
    temp = _mm512_xor_si512(vec1_temp, vec2_temp);
    temp = _mm512_and_si512(temp, mask);
    vec1_temp = _mm512_xor_si512(vec1_temp, temp);
    vec2_temp = _mm512_xor_si512(vec2_temp, temp);
    std::memcpy(vec1, &vec1_temp, 64);
    std::memcpy(vec2, &vec2_temp, 64);

    // on skylake, it's faster to run 256bit instructions two times
#else
    CXCHG_internal<32>(cond, vec1, vec2);
    CXCHG_internal<32>(cond, (char*)vec1 + 32, (char*)vec2 + 32);
#endif
    return;
  }
  if constexpr (sz >= 32) {
#if defined(__AVX512VL__)
    __m256d vec1_temp, vec2_temp;
    std::memcpy(&vec1_temp, vec1, 32);
    std::memcpy(&vec2_temp, vec2, 32);
    const __m256d& vec1_after_swap =
        _mm256_mask_blend_pd(blend_mask, vec1_temp, vec2_temp);
    const __m256d& vec2_after_swap =
        _mm256_mask_blend_pd(blend_mask, vec2_temp, vec1_temp);
    std::memcpy(vec1, &vec1_after_swap, 32);
    std::memcpy(vec2, &vec2_after_swap, 32);
#elif defined(__AVX2__)
    /* alternative implementation
    __m256i vec1_temp, vec2_temp;
    std::memcpy(&vec1_temp, vec1, 32);
    std::memcpy(&vec2_temp, vec2, 32);
    __m256i mask = _mm256_set1_epi8(-cond);
    const __m256i& vec1_after_swap =
        _mm256_blendv_epi8(vec1_temp, vec2_temp, mask);
    const __m256i& vec2_after_swap =
        _mm256_blendv_epi8(vec2_temp, vec1_temp, mask);
    std::memcpy(vec1, &vec1_after_swap, 32);
    std::memcpy(vec2, &vec2_after_swap, 32);
    */
    __m256i vec1_temp, vec2_temp;
    __m256i temp;
    std::memcpy(&vec1_temp, vec1, 32);
    std::memcpy(&vec2_temp, vec2, 32);
    __m256i mask = _mm256_set1_epi32(-cond);  // Create a mask based on cond
    temp = _mm256_xor_si256(vec1_temp, vec2_temp);
    temp = _mm256_and_si256(temp, mask);
    vec1_temp = _mm256_xor_si256(vec1_temp, temp);
    vec2_temp = _mm256_xor_si256(vec2_temp, temp);
    std::memcpy(vec1, &vec1_temp, 32);
    std::memcpy(vec2, &vec2_temp, 32);
#else
    CXCHG_internal<16>(cond, vec1, vec2);
    CXCHG_internal<16>(cond, (char*)vec1 + 16, (char*)vec2 + 16);
#endif
  }
  if constexpr (sz % 32 >= 16) {
    constexpr uint64_t offset = 4 * (sz / 32);
#if defined(__AVX512VL__)
    __m128d vec1_temp, vec2_temp;
    std::memcpy(&vec1_temp, (uint64_t*)vec1 + offset, 16);
    std::memcpy(&vec2_temp, (uint64_t*)vec2 + offset, 16);
    const __m128d& vec1_after_swap =
        _mm_mask_blend_pd(blend_mask, vec1_temp, vec2_temp);
    const __m128d& vec2_after_swap =
        _mm_mask_blend_pd(blend_mask, vec2_temp, vec1_temp);
    std::memcpy((uint64_t*)vec1 + offset, &vec1_after_swap, 16);
    std::memcpy((uint64_t*)vec2 + offset, &vec2_after_swap, 16);
#elif defined(__SSE2__)
    __m128i vec1_temp, vec2_temp;
    __m128i temp;
    std::memcpy(&vec1_temp, (uint64_t*)vec1 + offset, 16);
    std::memcpy(&vec2_temp, (uint64_t*)vec2 + offset, 16);
    __m128i mask = _mm_set1_epi16(-cond);  // Create a mask based on cond
    temp = _mm_xor_si128(vec1_temp, vec2_temp);
    temp = _mm_and_si128(temp, mask);
    vec1_temp = _mm_xor_si128(vec1_temp, temp);
    vec2_temp = _mm_xor_si128(vec2_temp, temp);
    std::memcpy((uint64_t*)vec1 + offset, &vec1_temp, 16);
    std::memcpy((uint64_t*)vec2 + offset, &vec2_temp, 16);
#else
    CXCHG_internal<8>(cond, (uint64_t*)vec1 + offset, (uint64_t*)vec2 + offset);
    CXCHG_internal<8>(cond, (char*)vec1 + 8 * offset + 8,
                      (char*)vec2 + 8 * offset + 8);
#endif
  }

  if constexpr (sz % 16 >= 8) {
    constexpr uint64_t offset = 2 * (sz / 16);
    uint64_t* curr1_64 = (uint64_t*)vec1 + offset;
    uint64_t* curr2_64 = (uint64_t*)vec2 + offset;
    CXCHG(cond, *curr1_64, *curr2_64);
  }
  if constexpr (sz % 8 >= 4) {
    constexpr uint64_t offset = 2 * (sz / 8);
    uint32_t* curr1_32 = (uint32_t*)vec1 + offset;
    uint32_t* curr2_32 = (uint32_t*)vec2 + offset;
    CXCHG(cond, *curr1_32, *curr2_32);
  }
  if constexpr (sz % 4 >= 2) {
    constexpr uint64_t offset = 2 * (sz / 4);
    uint16_t* curr1_16 = (uint16_t*)vec1 + offset;
    uint16_t* curr2_16 = (uint16_t*)vec2 + offset;
    CXCHG(cond, *curr1_16, *curr2_16);
  }
  if constexpr (sz % 2 >= 1) {
    constexpr uint64_t offset = 2 * (sz / 2);
    uint8_t* curr1_8 = (uint8_t*)vec1 + offset;
    uint8_t* curr2_8 = (uint8_t*)vec2 + offset;
    CXCHG(cond, *curr1_8, *curr2_8);
  }
}

template <const uint64_t sz>
INLINE void CMOV_internal(const bool cond, void* dest, const void* src) {
  static_assert(sz <= 64);
  if constexpr (sz == 64) {
#if defined(__AVX512VL__)

    const __mmask8 mask = (__mmask8)(!cond) - 1;
    // Create a mask based on cond
    __m512i srcVec = _mm512_loadu_si512(src);
    __m512i destVec =
        _mm512_mask_mov_epi64(_mm512_loadu_si512(dest), mask, srcVec);
    _mm512_storeu_si512(dest, destVec);

    // on skylake, it's faster to run 256bit instructions two times
#else
    CMOV_internal<32>(cond, dest, src);
    CMOV_internal<32>(cond, (char*)dest + 32, (char*)src + 32);
#endif
    return;
  }
  if constexpr (sz >= 32) {
#if defined(__AVX2__)
    const __m256i mask =
        _mm256_set1_epi64x(-(!!cond));  // Create a mask based on cond
    __m256i srcVec = _mm256_loadu_si256((__m256i*)src);
    __m256i destVec =
        _mm256_blendv_epi8(_mm256_loadu_si256((__m256i*)dest), srcVec, mask);
    _mm256_storeu_si256((__m256i*)dest, destVec);
#else
    CMOV_internal<16>(cond, dest, src);
    CMOV_internal<16>(cond, (char*)dest + 16, (char*)src + 16);
#endif
  }
  if constexpr (sz % 32 >= 16) {
    constexpr uint64_t offset = 4 * (sz / 32);
    uint64_t* src8 = (uint64_t*)src + offset;
    uint64_t* dest8 = (uint64_t*)dest + offset;
#if defined(__SSE2__)
    const __m128i mask =
        _mm_set1_epi64x(-(!!cond));  // Create a mask based on cond
    __m128i srcVec = _mm_loadu_si128((__m128i*)src8);
    __m128i destVec = _mm_loadu_si128((__m128i*)dest8);
    __m128i blended = _mm_or_si128(_mm_and_si128(mask, srcVec),
                                   _mm_andnot_si128(mask, destVec));
    _mm_storeu_si128((__m128i*)dest8, blended);
#else
    CMOV_internal<8>(cond, dest8, src8);
    CMOV_internal<8>(cond, (char*)dest8 + 8, (char*)src8 + 8);
#endif
  }

  if constexpr (sz % 16 >= 8) {
    constexpr uint64_t offset = 2 * (sz / 16);
    uint64_t* curr1_64 = (uint64_t*)dest + offset;
    uint64_t* curr2_64 = (uint64_t*)src + offset;
    CMOV(cond, *curr1_64, *curr2_64);
  }
  if constexpr (sz % 8 >= 4) {
    constexpr uint64_t offset = 2 * (sz / 8);
    uint32_t* curr1_32 = (uint32_t*)dest + offset;
    uint32_t* curr2_32 = (uint32_t*)src + offset;
    CMOV(cond, *curr1_32, *curr2_32);
  }
  if constexpr (sz % 4 >= 2) {
    constexpr uint64_t offset = 2 * (sz / 4);
    uint16_t* curr1_16 = (uint16_t*)dest + offset;
    uint16_t* curr2_16 = (uint16_t*)src + offset;
    CMOV(cond, *curr1_16, *curr2_16);
  }
  if constexpr (sz % 2 >= 1) {
    constexpr uint64_t offset = 2 * (sz / 2);
    uint8_t* curr1_8 = (uint8_t*)dest + offset;
    uint8_t* curr2_8 = (uint8_t*)src + offset;
    CMOV(cond, *curr1_8, *curr2_8);
  }
}

template <typename T>
INLINE void obliSwap(const bool mov, T& guy1, T& guy2) {
  // static_assert(sizeof(T)%8 == 0);
  __m512i* curr1 = (__m512i*)&guy1;
  __m512i* curr2 = (__m512i*)&guy2;
  for (uint64_t i = 0; i < sizeof(T) / 64; ++i) {
    CXCHG_internal<64>(mov, curr1, curr2);
    curr1++;
    curr2++;
  }
  constexpr uint64_t rem_size = sizeof(T) % 64;
  if constexpr (rem_size > 0) {
    CXCHG_internal<rem_size>(mov, curr1, curr2);
  }
}

// copy src to dst if mov is true, return mov
template <typename T>
INLINE bool obliMove(const bool mov, T& dest, const T& src) {
  __m512i* curr1 = (__m512i*)&dest;
  const __m512i* curr2 = (__m512i*)&src;
  for (uint64_t i = 0; i < sizeof(T) / 64; ++i) {
    CMOV_internal<64>(mov, curr1, curr2);
    curr1++;
    curr2++;
  }
  constexpr uint64_t rem_size = sizeof(T) % 64;
  if constexpr (rem_size > 0) {
    CMOV_internal<rem_size>(mov, curr1, curr2);
  }
  return mov;
}

// Some CMOV specializations:
//
template <>
INLINE void CMOV<uint64_t>(const uint64_t& cond, uint64_t& val1,
                           const uint64_t& val2) {
  CMOV8(cond, val1, val2);
}

template <>
INLINE void CMOV<uint32_t>(const uint64_t& cond, uint32_t& val1,
                           const uint32_t& val2) {
  CMOV4(cond, val1, val2);
}

template <>
INLINE void CMOV<uint16_t>(const uint64_t& cond, uint16_t& val1,
                           const uint16_t& val2) {
  CMOV2(cond, val1, val2);
}

template <>
INLINE void CMOV<uint8_t>(const uint64_t& cond, uint8_t& val1,
                          const uint8_t& val2) {
  CMOV1(cond, val1, val2);
}

template <>
INLINE void CMOV<bool>(const uint64_t& cond, bool& val1, const bool& val2) {
  CMOV_BOOL(cond, val1, val2);
}

template <>
INLINE void CMOV<int>(const uint64_t& cond, int& val1, const int& val2) {
  CMOV4(cond, (uint32_t&)val1, val2);
}

template <>
INLINE void CMOV<short>(const uint64_t& cond, short& val1, const short& val2) {
  CMOV2(cond, (uint16_t&)val1, val2);
}

template <>
INLINE void CMOV<int8_t>(const uint64_t& cond, int8_t& val1,
                         const int8_t& val2) {
  CMOV1(cond, (uint8_t&)val1, val2);
}

// Other cmov specializations are inside the respective headers for the types.
// Make sure to always include this file first so that overloads get resolved,
// correctly.
//
// Aditionally, TSET and CXCHG need to be overloaded if CMOV is overloaded
// for a given type also:
//

#define OVERLOAD_TSET_CXCHG(TYPE, ...)                                     \
                                                                           \
  template <__VA_ARGS__>                                                   \
  INLINE void TSET(bool selector, TYPE& A, const TYPE& B, const TYPE& C) { \
    CMOV(selector, A, C);                                                  \
    CMOV(!selector, A, B);                                                 \
  }                                                                        \
                                                                           \
  template <__VA_ARGS__>                                                   \
  INLINE void CTSET(bool condition, bool selector, TYPE& A, const TYPE& B, \
                    const TYPE& C) {                                       \
    CMOV(condition* selector, A, C);                                       \
    CMOV(condition * !selector, A, B);                                     \
  }                                                                        \
                                                                           \
  template <__VA_ARGS__>                                                   \
  INLINE void CXCHG(const uint64_t& cond, TYPE& A, TYPE& B) {              \
    if constexpr (sizeof(TYPE(A)) == 8) {                                  \
      CSWAP8(cond, (uint64_t&)A, (uint64_t&)B);                            \
      return;                                                              \
    }                                                                      \
    const TYPE C = A;                                                      \
    CMOV(cond, A, B);                                                      \
    CMOV(cond, B, C);                                                      \
  }

// Overload for CMOV of pair
template <typename A, typename B>
INLINE void CMOV(const uint64_t& cond, std::pair<A, B>& val1,
                 const std::pair<A, B>& val2) {
  CMOV(cond, val1.first, val2.first);
  CMOV(cond, val1.second, val2.second);
}

OVERLOAD_TSET_CXCHG(std::pair<X COMMA Y>, typename X, typename Y);

#ifdef __AVX2__
#define m256i __m256i
#define mm256_set1_epi32 _mm256_set1_epi32
#define mm256_contain_le_zero(vec) \
  _mm256_movemask_epi8(_mm256_cmpgt_epi32(_mm256_set1_epi32(1), vec))
INLINE uint32_t mm256_extract_epi32_var_indx(const m256i vec,
                                             const unsigned int i) {
  __m128i indx = _mm_cvtsi32_si128(i);
  __m256i val = _mm256_permutevar8x32_epi32(vec, _mm256_castsi128_si256(indx));
  return _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}

// decrement vec at index i, if i >= 8, return vec
INLINE m256i mm256_decrement_epi32_var_indx(const m256i vec,
                                            const unsigned int i) {
  static const __m256i mask = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
  __m256i cmp = _mm256_set1_epi32(
      i);  // create 256-bit register with i in all 32-bit lanes
  __m256i cmp_result = _mm256_cmpeq_epi32(
      mask, cmp);  // will set the matching 32 bits as 111..11
  return _mm256_add_epi32(vec, cmp_result);  // decrement
}
#else
struct m256i {
  int32_t data[8];
};
static INLINE int32_t mm256_extract_epi32_var_indx(const m256i vec,
                                                   const unsigned int i) {
  int32_t ans;
  for (unsigned int j = 0; j < 8; ++j) {
    CMOV(i == j, ans, vec.data[j]);
  }
  return ans;
}
static INLINE m256i mm256_decrement_epi32_var_indx(m256i vec,
                                                   const unsigned int i) {
  for (unsigned int j = 0; j < 8; ++j) {
    CMOV(i == j, vec.data[j], vec.data[j] - 1);
  }
  return vec;
}
static INLINE bool mm256_contain_le_zero(m256i vec) {
  bool res = false;
  for (unsigned int j = 0; j < 8; ++j) {
    res |= vec.data[j] <= 0;
  }
  return res;
}
static INLINE m256i mm256_set1_epi32(int32_t z) {
  m256i res;
  for (unsigned int j = 0; j < 8; ++j) {
    res.data[j] = z;
  }
  return res;
}

#endif