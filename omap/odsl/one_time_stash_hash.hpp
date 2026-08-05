#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "algorithm/bitonic.hpp"
#include "algorithm/or_compact_shuffle.hpp"
#include "common/dummy.hpp"
#include "common/encutils.hpp"
#include "common/mov_intrinsics.hpp"
#include "common/tracing/perf.hpp"

namespace ODSL {

namespace detail {

template <typename T>
INLINE void SecureWipe(T& value) {
  volatile uint8_t* bytes = reinterpret_cast<volatile uint8_t*>(&value);
  for (size_t i = 0; i < sizeof(T); ++i) {
    bytes[i] = 0;
  }
}

template <typename T>
INLINE void SecureWipe(std::vector<T>& values) {
  volatile uint8_t* bytes = reinterpret_cast<volatile uint8_t*>(values.data());
  for (size_t i = 0; i < values.size() * sizeof(T); ++i) {
    bytes[i] = 0;
  }
}

}  // namespace detail

/**
 * A separately keyed AES PRF for one-time stash tables.
 *
 * The 96-bit CTR nonce is the canonical big-endian encoding of a 32-bit
 * domain followed by a 64-bit value. Encrypting a zero block therefore gives
 * AES(key, domain || value || 0). The context and key are local to one helper;
 * neither the external-memory encryption key nor a long-lived OMap salt is
 * reused.
 */
template <typename Uid>
class OneTimeStashAesHash {
 public:
  using Key = std::array<uint8_t, 32>;

 private:
#if defined(__AES__)
  using Context = br_aes_x86ni_ctr_keys;
#else
  using Context = br_aes_ct_ctr_keys;
#endif

  Key key_{};
  Context context_{};

#if defined(__VAES__) && defined(__AVX512F__)
  struct alignas(64) WideRoundKey {
    __m512i value;
  };
  std::array<WideRoundKey, 15> wideRoundKeys_{};
#endif

  INLINE static void EncodeInput(uint64_t value, uint32_t domain,
                                 uint8_t* block) {
    for (size_t i = 0; i < 4; ++i) {
      block[i] = static_cast<uint8_t>(domain >> (8 * (3 - i)));
    }
    for (size_t i = 0; i < 8; ++i) {
      block[4 + i] = static_cast<uint8_t>(value >> (8 * (7 - i)));
    }
    std::memset(block + 12, 0, 4);
  }

  INLINE static uint64_t DecodeResult(const uint8_t* block) {
    uint64_t result = 0;
    for (size_t i = 0; i < sizeof(result); ++i) {
      result = (result << 8) | block[i];
    }
    return result;
  }

#if defined(__AES__)
  INLINE void EncryptFourBlocks(uint8_t* blocks) const {
#if defined(__VAES__) && defined(__AVX512F__)
    __m512i state = _mm512_load_si512(blocks);
    state = _mm512_xor_si512(state, wideRoundKeys_[0].value);
#pragma GCC unroll 14
    for (size_t round = 1; round < 14; ++round) {
      state = _mm512_aesenc_epi128(state, wideRoundKeys_[round].value);
    }
    state = _mm512_aesenclast_epi128(state, wideRoundKeys_[14].value);
    _mm512_store_si512(blocks, state);
#else
    __m128i states[4];
    const __m128i firstKey = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(context_.skey.skni));
    for (size_t lane = 0; lane < 4; ++lane) {
      states[lane] = _mm_xor_si128(
          _mm_load_si128(reinterpret_cast<const __m128i*>(blocks + 16 * lane)),
          firstKey);
    }
    for (size_t round = 1; round < 14; ++round) {
      const __m128i roundKey = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(context_.skey.skni + 16 * round));
      for (auto& state : states) {
        state = _mm_aesenc_si128(state, roundKey);
      }
    }
    const __m128i lastKey = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(context_.skey.skni + 16 * 14));
    for (size_t lane = 0; lane < 4; ++lane) {
      states[lane] = _mm_aesenclast_si128(states[lane], lastKey);
      _mm_store_si128(reinterpret_cast<__m128i*>(blocks + 16 * lane),
                      states[lane]);
    }
#endif
  }
#endif

  template <typename GetValue>
  void HashValuesToBuckets(size_t count, uint32_t domain, uint32_t modulus,
                           uint32_t* buckets, GetValue&& getValue) const {
    Assert(modulus != 0);
#if defined(__AES__)
    alignas(64) uint8_t blocks[64];
    for (size_t offset = 0; offset < count; offset += 4) {
      const size_t lanes = std::min<size_t>(4, count - offset);
      std::memset(blocks, 0, sizeof(blocks));
      for (size_t lane = 0; lane < lanes; ++lane) {
        EncodeInput(getValue(offset + lane), domain, blocks + 16 * lane);
      }
      EncryptFourBlocks(blocks);
      for (size_t lane = 0; lane < lanes; ++lane) {
        buckets[offset + lane] = static_cast<uint32_t>(
            DecodeResult(blocks + 16 * lane) % modulus);
      }
    }
    detail::SecureWipe(blocks);
#else
    for (size_t i = 0; i < count; ++i) {
      buckets[i] = static_cast<uint32_t>(HashValue(getValue(i), domain) %
                                         modulus);
    }
#endif
  }

  INLINE uint64_t HashValue(uint64_t value, uint32_t domain) const {
    uint8_t nonce[12];
    uint8_t input[16];
    EncodeInput(value, domain, input);
    std::memcpy(nonce, input, sizeof(nonce));

    uint8_t tag[16] = {0};
#if defined(__AES__)
    br_aes_x86ni_ctr_run(&context_, nonce, 0, tag, sizeof(tag));
#else
    br_aes_ct_ctr_run(&context_, nonce, 0, tag, sizeof(tag));
#endif
    const uint64_t result = DecodeResult(tag);
    detail::SecureWipe(tag);
    detail::SecureWipe(nonce);
    detail::SecureWipe(input);
    return result;
  }

 public:
  OneTimeStashAesHash() = default;
  ~OneTimeStashAesHash() { Wipe(); }

  void GenerateKey(Key& key) { read_rand(key.data(), key.size()); }

  void SetKey(const Key& key) {
    key_ = key;
#if defined(__AES__)
    br_aes_x86ni_ctr_init(&context_, key_.data(), key_.size());
#if defined(__VAES__) && defined(__AVX512F__)
    for (size_t round = 0; round < wideRoundKeys_.size(); ++round) {
      const __m128i roundKey = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(context_.skey.skni + 16 * round));
      wideRoundKeys_[round].value = _mm512_broadcast_i32x4(roundKey);
    }
#endif
#else
    br_aes_ct_ctr_init(&context_, key_.data(), key_.size());
#endif
  }

  INLINE uint64_t HashUid(const Uid& uid, uint32_t domain) const {
    static_assert(std::is_integral_v<Uid> && std::is_unsigned_v<Uid>,
                  "one-time stash hashing requires an unsigned integer UID");
    static_assert(sizeof(Uid) <= sizeof(uint64_t),
                  "one-time stash hashing supports UIDs up to 64 bits");
    return HashValue(static_cast<uint64_t>(uid), domain);
  }

  INLINE uint64_t HashIndex(uint64_t requestIndex, uint32_t domain) const {
    return HashValue(requestIndex, domain);
  }

  void HashUidBuckets(const Uid* uids, size_t count, uint32_t domain,
                      uint32_t modulus, uint32_t* buckets) const {
    HashValuesToBuckets(count, domain, modulus, buckets,
                        [&](size_t i) {
                          return static_cast<uint64_t>(uids[i]);
                        });
  }

  void HashIndexBuckets(size_t count, uint32_t domain, uint32_t modulus,
                        uint32_t* buckets) const {
    HashValuesToBuckets(count, domain, modulus, buckets,
                        [](size_t i) { return static_cast<uint64_t>(i); });
  }

  void Wipe() {
#if defined(__VAES__) && defined(__AVX512F__)
    detail::SecureWipe(wideRoundKeys_);
#endif
    detail::SecureWipe(context_);
    detail::SecureWipe(key_);
  }
};

/**
 * A freshly keyed, one-choice bucketed representation of a flat ORAM stash.
 *
 * Build leaves the input stash untouched until a candidate succeeds. Once it
 * succeeds, this object is authoritative until Restore is called. Destruction
 * also restores an active table, so ordinary C++ exceptions cannot expose a
 * stale flat stash to later ORAM operations.
 *
 * HashProvider is injectable for deterministic collision tests. It must
 * provide Key, GenerateKey, SetKey, HashUid, HashIndex, and Wipe.
 */
template <typename Block, typename Uid, size_t S, size_t B, size_t M, size_t G,
          size_t Candidates = 1,
          typename HashProvider = OneTimeStashAesHash<Uid>>
class OneTimeStashHash {
 public:
  static constexpr size_t MainSlots = M * B;
  static constexpr size_t SlotCount = MainSlots + G;
  static constexpr uint32_t MainDomain = 0x53544d31;       // "STM1"
  static constexpr uint32_t DuplicateDomain = 0x53544431;  // "STD1"

  using Key = typename HashProvider::Key;

 private:
  static_assert(S > 0, "stash size must be positive");
  static_assert(B > 0, "main bucket size must be positive");
  static_assert(M > 0, "the table must have a main bucket");
  static_assert(G > 0, "the table must have an overflow bucket");
  static_assert(Candidates > 0, "at least one candidate is required");
  static_assert(S <= std::numeric_limits<uint32_t>::max());
  static_assert(B <= std::numeric_limits<uint32_t>::max());
  static_assert(M <= std::numeric_limits<uint32_t>::max());
  static_assert(SlotCount < std::numeric_limits<uint32_t>::max());
  static_assert(SlotCount >= S,
                "the one-time table must have at least S total slots");
  static_assert(std::is_trivially_copyable_v<Block>);
  static_assert(std::is_trivially_copyable_v<Key>);

  static constexpr size_t CounterLanes = ((M + 15) / 16) * 16;

  HashProvider hasher_;
  Key selectedKey_{};
  Block* flatStash_ = nullptr;
  bool active_ = false;
  bool workspaceWiped_ = false;

  // Reused, fixed-size workspaces. Keeping marks_ allocated makes the
  // destructor's restoration path non-allocating.
  std::array<uint32_t, CounterLanes> counters_{};
  std::array<uint32_t, CounterLanes> counterLaneIds_{};
  std::array<uint32_t, S> destinations_{};
  std::array<Uid, S> buildUids_{};
  std::array<uint32_t, S> buildBuckets_{};
  std::vector<Block> records_;
  std::vector<Block> candidate_;
  std::vector<Block> table_;
  std::vector<uint32_t> marks_;
  std::vector<uint32_t> queryRealBuckets_;
  std::vector<uint32_t> queryDuplicateBuckets_;

  INLINE static Block DummyBlock() { return Block(); }

  static void FillUidBuckets(const HashProvider& hasher, const Uid* uids,
                             size_t count, uint32_t domain,
                             uint32_t* buckets) {
    if constexpr (requires {
                    hasher.HashUidBuckets(uids, count, domain,
                                          static_cast<uint32_t>(M), buckets);
                  }) {
      hasher.HashUidBuckets(uids, count, domain, static_cast<uint32_t>(M),
                            buckets);
    } else {
      for (size_t i = 0; i < count; ++i) {
        buckets[i] = static_cast<uint32_t>(hasher.HashUid(uids[i], domain) % M);
      }
    }
  }

  static void FillIndexBuckets(const HashProvider& hasher, size_t count,
                               uint32_t domain, uint32_t* buckets) {
    if constexpr (requires {
                    hasher.HashIndexBuckets(count, domain,
                                            static_cast<uint32_t>(M), buckets);
                  }) {
      hasher.HashIndexBuckets(count, domain, static_cast<uint32_t>(M),
                              buckets);
    } else {
      for (size_t i = 0; i < count; ++i) {
        buckets[i] =
            static_cast<uint32_t>(hasher.HashIndex(i, domain) % M);
      }
    }
  }

  /** Increment one secret counter while scanning every public counter lane. */
  INLINE uint32_t IncrementAndGetRank(uint32_t bucket, bool isReal) {
    uint32_t rank = 0;

#if defined(__AVX512F__)
    const __m512i target = _mm512_set1_epi32(static_cast<int>(bucket));
    const __m512i one = _mm512_set1_epi32(1);
    for (size_t offset = 0; offset < CounterLanes; offset += 16) {
      const __m512i lanes = _mm512_loadu_si512(counterLaneIds_.data() + offset);
      const __m512i oldCounters = _mm512_loadu_si512(counters_.data() + offset);
      __mmask16 matches = _mm512_cmpeq_epi32_mask(lanes, target);
      matches &= static_cast<__mmask16>(-static_cast<int>(isReal));
      const __m512i newCounters =
          _mm512_mask_add_epi32(oldCounters, matches, oldCounters, one);
      _mm512_storeu_si512(counters_.data() + offset, newCounters);

      alignas(64) uint32_t oldLanes[16];
      _mm512_store_si512(oldLanes, oldCounters);
      for (size_t lane = 0; lane < 16; ++lane) {
        const bool match = isReal & (bucket == offset + lane);
        obliMove(match, rank, oldLanes[lane] + 1);
      }
    }
#elif defined(__AVX2__)
    const __m256i target = _mm256_set1_epi32(static_cast<int>(bucket));
    const __m256i one = _mm256_set1_epi32(1);
    const __m256i realMask = _mm256_set1_epi32(-static_cast<int>(isReal));
    for (size_t offset = 0; offset < CounterLanes; offset += 8) {
      const __m256i lanes = _mm256_loadu_si256(
          reinterpret_cast<const __m256i*>(counterLaneIds_.data() + offset));
      const __m256i oldCounters = _mm256_loadu_si256(
          reinterpret_cast<const __m256i*>(counters_.data() + offset));
      const __m256i matches =
          _mm256_and_si256(_mm256_cmpeq_epi32(lanes, target), realMask);
      const __m256i newCounters =
          _mm256_add_epi32(oldCounters, _mm256_and_si256(matches, one));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(counters_.data() + offset),
                          newCounters);

      alignas(32) uint32_t oldLanes[8];
      _mm256_store_si256(reinterpret_cast<__m256i*>(oldLanes), oldCounters);
      for (size_t lane = 0; lane < 8; ++lane) {
        const bool match = isReal & (bucket == offset + lane);
        obliMove(match, rank, oldLanes[lane] + 1);
      }
    }
#else
    for (size_t lane = 0; lane < CounterLanes; ++lane) {
      const bool match = isReal & (bucket == lane);
      obliMove(match, rank, counters_[lane] + 1);
      counters_[lane] += match;
    }
#endif

    return rank;
  }

  bool BuildCandidate(const Block* flatStash, HashProvider& candidateHasher) {
    std::fill(counters_.begin(), counters_.end(), 0);
    std::fill(candidate_.begin(), candidate_.end(), DummyBlock());

    for (size_t i = 0; i < S; ++i) {
      buildUids_[i] = flatStash[i].uid;
    }
    FillUidBuckets(candidateHasher, buildUids_.data(), S, MainDomain,
                   buildBuckets_.data());

    uint32_t residualCount = 0;
    bool failed = false;
    for (size_t i = 0; i < S; ++i) {
      records_[i] = flatStash[i];
      const bool isReal = !records_[i].IsDummy();
      const uint32_t bucket = buildBuckets_[i];
      const uint32_t rank = IncrementAndGetRank(bucket, isReal);
      const bool isResidual = isReal & (rank > B);
      const uint32_t residualRank = residualCount;
      residualCount += isResidual;

      uint32_t destination = static_cast<uint32_t>(SlotCount);
      const uint32_t mainDestination =
          bucket * static_cast<uint32_t>(B) + rank - 1;
      const uint32_t residualDestination =
          static_cast<uint32_t>(MainSlots) + residualRank;
      obliMove(isReal & !isResidual, destination, mainDestination);
      obliMove(isResidual & (residualRank < G), destination,
               residualDestination);
      destinations_[i] = destination;

      const bool overflow = isResidual & (residualRank >= G);
      failed |= overflow;
      const Uid dummy = DUMMY<Uid>();
      obliMove(overflow, records_[i].uid, dummy);
    }

    Algorithm::BitonicSortSepPayload(destinations_.begin(), destinations_.end(),
                                     records_.begin());
    for (size_t i = 0; i < S; ++i) {
      candidate_[i] = records_[i];
    }

    uint32_t prefix = 0;
    marks_[0] = 0;
    size_t slot = 0;
    for (size_t bucket = 0; bucket < M; ++bucket) {
      for (size_t rank = 0; rank < B; ++rank, ++slot) {
        prefix += counters_[bucket] > rank;
        marks_[slot + 1] = prefix;
      }
    }
    for (size_t rank = 0; rank < G; ++rank, ++slot) {
      prefix += residualCount > rank;
      marks_[slot + 1] = prefix;
    }
    Algorithm::OrDistributeSeparateMark(candidate_.begin(), candidate_.end(),
                                        marks_.begin());
    return !failed;
  }

  void WipeWorkspace() {
    detail::SecureWipe(selectedKey_);
    detail::SecureWipe(counters_);
    detail::SecureWipe(destinations_);
    detail::SecureWipe(buildUids_);
    detail::SecureWipe(buildBuckets_);
    detail::SecureWipe(records_);
    detail::SecureWipe(candidate_);
    detail::SecureWipe(table_);
    detail::SecureWipe(marks_);
    detail::SecureWipe(queryRealBuckets_);
    detail::SecureWipe(queryDuplicateBuckets_);
    hasher_.Wipe();
    workspaceWiped_ = true;
  }

 public:
  explicit OneTimeStashHash(HashProvider hasher = HashProvider())
      : hasher_(std::move(hasher)),
        records_(S),
        candidate_(SlotCount),
        table_(SlotCount),
        marks_(SlotCount + 1) {
    for (size_t i = 0; i < CounterLanes; ++i) {
      counterLaneIds_[i] = static_cast<uint32_t>(i);
    }
  }

  OneTimeStashHash(const OneTimeStashHash&) = delete;
  OneTimeStashHash& operator=(const OneTimeStashHash&) = delete;

  ~OneTimeStashHash() noexcept {
    if (active_) {
      Restore();
    } else if (!workspaceWiped_) {
      WipeWorkspace();
    }
    detail::SecureWipe(counterLaneIds_);
  }

  /** Build all candidates and select the first successful one obliviously. */
  bool Build(Block* flatStash) {
    Assert(!active_);
    workspaceWiped_ = false;
    flatStash_ = flatStash;
    std::fill(table_.begin(), table_.end(), DummyBlock());
    detail::SecureWipe(selectedKey_);

    bool selected = false;
    for (size_t candidateIndex = 0; candidateIndex < Candidates;
         ++candidateIndex) {
      Key candidateKey{};
      hasher_.GenerateKey(candidateKey);
      HashProvider candidateHasher = hasher_;
      candidateHasher.SetKey(candidateKey);
      const bool success = BuildCandidate(flatStash, candidateHasher);
      const bool select = success & !selected;
      for (size_t slot = 0; slot < SlotCount; ++slot) {
        obliMove(select, table_[slot], candidate_[slot]);
      }
      uint8_t* selectedKeyBytes = reinterpret_cast<uint8_t*>(&selectedKey_);
      const uint8_t* candidateKeyBytes =
          reinterpret_cast<const uint8_t*>(&candidateKey);
      for (size_t i = 0; i < sizeof(Key); ++i) {
        obliMove(select, selectedKeyBytes[i], candidateKeyBytes[i]);
      }
      selected |= success;
      PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_CANDIDATE_OVERFLOWS,
                           !success);
      candidateHasher.Wipe();
      detail::SecureWipe(candidateKey);
    }

    PERFCTR_INCREMENT(CIRCUITORAM_STASH_HASH_BUILDS);
    PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_FALLBACKS, !selected);
    active_ = selected;
    if (selected) {
      hasher_.SetKey(selectedKey_);
    } else {
      flatStash_ = nullptr;
      WipeWorkspace();
    }
    return selected;
  }

  /** Precompute the real and duplicate-probe buckets for one request batch. */
  void PrepareBatch(const Uid* uids, size_t count) {
    Assert(active_);
    Assert(count == 0 || uids != nullptr);
    queryRealBuckets_.resize(count);
    queryDuplicateBuckets_.resize(count);
    FillUidBuckets(hasher_, uids, count, MainDomain,
                   queryRealBuckets_.data());
    FillIndexBuckets(hasher_, count, DuplicateDomain,
                     queryDuplicateBuckets_.data());
  }

  /** Probe exactly one prepared main bucket and the shared overflow bucket. */
  template <typename Data>
  bool ReadAndRemovePrepared(const Uid& uid, size_t requestIndex,
                             bool firstOccurrence, Data& out) {
    Assert(active_);
    Assert(requestIndex < queryRealBuckets_.size());
    uint32_t bucket = queryDuplicateBuckets_[requestIndex];
    obliMove(firstOccurrence, bucket, queryRealBuckets_[requestIndex]);

    bool found = false;
    const size_t mainOffset = static_cast<size_t>(bucket) * B;
    for (size_t i = 0; i < B; ++i) {
      found |= table_[mainOffset + i].ReadAndRemove(uid, out);
    }
    for (size_t i = 0; i < G; ++i) {
      found |= table_[MainSlots + i].ReadAndRemove(uid, out);
    }
    PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_MAIN_ENTRIES_PROBED, B);
    PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_OVERFLOW_ENTRIES_PROBED, G);
    return found;
  }

  /** Scalar compatibility path used by focused tests and uncommon callers. */
  template <typename Data>
  bool ReadAndRemove(const Uid& uid, uint64_t requestIndex,
                     bool firstOccurrence, Data& out) {
    Assert(active_);
    const uint32_t realBucket =
        static_cast<uint32_t>(hasher_.HashUid(uid, MainDomain) % M);
    const uint32_t duplicateBucket = static_cast<uint32_t>(
        hasher_.HashIndex(requestIndex, DuplicateDomain) % M);
    uint32_t bucket = duplicateBucket;
    obliMove(firstOccurrence, bucket, realBucket);

    bool found = false;
    const size_t mainOffset = static_cast<size_t>(bucket) * B;
    for (size_t i = 0; i < B; ++i) {
      found |= table_[mainOffset + i].ReadAndRemove(uid, out);
    }
    for (size_t i = 0; i < G; ++i) {
      found |= table_[MainSlots + i].ReadAndRemove(uid, out);
    }
    PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_MAIN_ENTRIES_PROBED, B);
    PERFCTR_INCREMENT_BY(CIRCUITORAM_STASH_HASH_OVERFLOW_ENTRIES_PROBED, G);
    return found;
  }

  /** Obliviously compact the table and restore the ordinary S-slot stash. */
  void Restore() noexcept {
    if (!active_) {
      return;
    }

    uint32_t prefix = 0;
    marks_[0] = 0;
    for (size_t i = 0; i < SlotCount; ++i) {
      prefix += !table_[i].IsDummy();
      marks_[i + 1] = prefix;
    }
    Algorithm::OrCompactSeparateMark(table_.begin(), table_.end(),
                                     marks_.begin());
    for (size_t i = 0; i < S; ++i) {
      flatStash_[i] = table_[i];
    }

    active_ = false;
    flatStash_ = nullptr;
    PERFCTR_INCREMENT(CIRCUITORAM_STASH_HASH_RESTORES);
    WipeWorkspace();
  }

  bool IsActive() const { return active_; }
};

}  // namespace ODSL
