#include "odsl/one_time_stash_hash.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "odsl/circuit_oram.hpp"

namespace {

template <typename Uid>
class DeterministicStashHash {
 public:
  using Key = uint64_t;

 private:
  Key nextKey_ = 0;
  Key key_ = 0;

 public:
  explicit DeterministicStashHash(Key firstKey = 0) : nextKey_(firstKey) {}

  void GenerateKey(Key& key) { key = nextKey_++; }
  void SetKey(const Key& key) { key_ = key; }

  uint64_t HashUid(const Uid& uid, uint32_t) const {
    // Key zero deliberately puts everything in one bucket. Later keys spread
    // consecutive test UIDs over consecutive buckets.
    return key_ == 0 ? 0 : static_cast<uint64_t>(uid);
  }

  uint64_t HashIndex(uint64_t requestIndex, uint32_t) const {
    return requestIndex * 0x9e3779b97f4a7c15ULL + key_;
  }

  void Wipe() { key_ = 0; }
};

using TestBlock = ODSL::Block<uint64_t, uint32_t, uint32_t>;

TestBlock MakeBlock(uint32_t uid) {
  return TestBlock(1000 + uid, 2000 + uid, uid);
}

TEST(OneTimeStashAesHash, BatchedBucketsMatchScalarPrf) {
  ODSL::OneTimeStashAesHash<uint64_t> hash;
  ODSL::OneTimeStashAesHash<uint64_t>::Key key;
  for (size_t i = 0; i < key.size(); ++i) {
    key[i] = static_cast<uint8_t>(3 * i + 1);
  }
  hash.SetKey(key);

  constexpr size_t count = 17;
  constexpr uint32_t modulus = 80;
  std::array<uint64_t, count> uids;
  std::array<uint32_t, count> uidBuckets;
  std::array<uint32_t, count> indexBuckets;
  for (size_t i = 0; i < count; ++i) {
    uids[i] = (uint64_t{1} << 48) + i * UINT64_C(0x9e3779b97f4a7c1);
  }

  hash.HashUidBuckets(uids.data(), uids.size(), 0x53544d31, modulus,
                      uidBuckets.data());
  hash.HashIndexBuckets(count, 0x53544431, modulus, indexBuckets.data());
  for (size_t i = 0; i < count; ++i) {
    EXPECT_EQ(uidBuckets[i], hash.HashUid(uids[i], 0x53544d31) % modulus);
    EXPECT_EQ(indexBuckets[i], hash.HashIndex(i, 0x53544431) % modulus);
  }
}

TEST(OneTimeStashHash, SharedOverflowBoundaryAndFailureFallback) {
  using Hash = ODSL::OneTimeStashHash<TestBlock, uint32_t, 9, 4, 8, 4, 1,
                                      DeterministicStashHash<uint32_t>>;

  std::array<TestBlock, 9> stash;
  for (size_t i = 0; i < 8; ++i) {
    stash[i] = MakeBlock(static_cast<uint32_t>(i + 1));
  }

  Hash atBoundary(DeterministicStashHash<uint32_t>(0));
  ASSERT_TRUE(atBoundary.Build(stash.data()));
  for (uint32_t uid = 1; uid <= 8; ++uid) {
    uint64_t value = 0;
    EXPECT_TRUE(atBoundary.ReadAndRemove(uid, uid - 1, true, value));
    EXPECT_EQ(value, 1000 + uid);
  }
  atBoundary.Restore();
  EXPECT_TRUE(
      std::all_of(stash.begin(), stash.end(),
                  [](const TestBlock& block) { return block.IsDummy(); }));

  for (size_t i = 0; i < stash.size(); ++i) {
    stash[i] = MakeBlock(static_cast<uint32_t>(i + 1));
  }
  const auto original = stash;
  Hash overflow(DeterministicStashHash<uint32_t>(0));
  EXPECT_FALSE(overflow.Build(stash.data()));
  for (size_t i = 0; i < stash.size(); ++i) {
    EXPECT_EQ(stash[i].uid, original[i].uid);
    EXPECT_EQ(stash[i].position, original[i].position);
    EXPECT_EQ(stash[i].data, original[i].data);
  }
}

TEST(OneTimeStashHash, SelectsLaterCandidateAndRestoresSurvivors) {
  using Hash = ODSL::OneTimeStashHash<TestBlock, uint32_t, 9, 4, 8, 4, 2,
                                      DeterministicStashHash<uint32_t>>;

  std::array<TestBlock, 9> stash;
  for (size_t i = 0; i < stash.size(); ++i) {
    stash[i] = MakeBlock(static_cast<uint32_t>(i + 1));
  }

  Hash hash(DeterministicStashHash<uint32_t>(0));
  ASSERT_TRUE(hash.Build(stash.data()));  // candidate 0 fails; candidate 1 fits

  uint64_t value = 0;
  EXPECT_TRUE(hash.ReadAndRemove(2, 0, true, value));
  EXPECT_EQ(value, 1002);
  EXPECT_FALSE(hash.ReadAndRemove(2, 1, false, value));
  EXPECT_TRUE(hash.ReadAndRemove(7, 2, true, value));
  EXPECT_EQ(value, 1007);
  EXPECT_FALSE(hash.ReadAndRemove(99, 3, true, value));
  hash.Restore();

  std::unordered_map<uint32_t, TestBlock> survivors;
  for (const auto& block : stash) {
    if (!block.IsDummy()) {
      survivors.emplace(block.uid, block);
    }
  }
  ASSERT_EQ(survivors.size(), 7UL);
  for (uint32_t uid = 1; uid <= 9; ++uid) {
    if (uid == 2 || uid == 7) {
      EXPECT_FALSE(survivors.contains(uid));
    } else {
      ASSERT_TRUE(survivors.contains(uid));
      EXPECT_EQ(survivors.at(uid).data, 1000 + uid);
      EXPECT_EQ(survivors.at(uid).position, 2000 + uid);
    }
  }
}

TEST(OneTimeStashHash, ResidualsAreSharedAcrossMainBuckets) {
  using Hash = ODSL::OneTimeStashHash<TestBlock, uint32_t, 14, 4, 8, 4, 1,
                                      DeterministicStashHash<uint32_t>>;

  std::array<TestBlock, 14> stash;
  for (uint32_t i = 0; i < 6; ++i) {
    stash[i] = MakeBlock(8 * (i + 1));
    stash[6 + i] = MakeBlock(8 * (i + 1) + 1);
  }
  Hash atBoundary(DeterministicStashHash<uint32_t>(1));
  ASSERT_TRUE(atBoundary.Build(stash.data()));
  atBoundary.Restore();

  stash[12] = MakeBlock(8 * 7);
  Hash overflow(DeterministicStashHash<uint32_t>(1));
  EXPECT_FALSE(overflow.Build(stash.data()));
}

TEST(OneTimeStashHash, Supports64BitUidAndPosition) {
  using Block64 = ODSL::Block<std::array<uint64_t, 4>, uint64_t, uint64_t>;
  using Hash = ODSL::OneTimeStashHash<Block64, uint64_t, 5, 4, 8, 4, 1,
                                      DeterministicStashHash<uint64_t>>;

  std::array<Block64, 5> stash;
  const uint64_t uid = (uint64_t{1} << 40) + 17;
  const std::array<uint64_t, 4> data = {1, 2, 3, 4};
  stash[3] = Block64(data, (uint64_t{1} << 48) + 9, uid);

  Hash hash(DeterministicStashHash<uint64_t>(1));
  ASSERT_TRUE(hash.Build(stash.data()));
  std::array<uint64_t, 4> out{};
  EXPECT_TRUE(hash.ReadAndRemove(uid, 0, true, out));
  EXPECT_EQ(out, data);
  hash.Restore();
  EXPECT_TRUE(std::all_of(stash.begin(), stash.end(), [](const Block64& block) {
    return block.IsDummy();
  }));
}

TEST(OneTimeStashHash, DestructorRestoresAuthoritativeTable) {
  using Hash = ODSL::OneTimeStashHash<TestBlock, uint32_t, 9, 4, 8, 4, 1,
                                      DeterministicStashHash<uint32_t>>;

  std::array<TestBlock, 9> stash;
  for (size_t i = 0; i < stash.size(); ++i) {
    stash[i] = MakeBlock(static_cast<uint32_t>(i + 1));
  }
  {
    Hash hash(DeterministicStashHash<uint32_t>(1));
    ASSERT_TRUE(hash.Build(stash.data()));
    uint64_t value = 0;
    ASSERT_TRUE(hash.ReadAndRemove(4, 0, true, value));
  }

  size_t realCount = 0;
  for (const auto& block : stash) {
    realCount += !block.IsDummy();
    EXPECT_TRUE(block.IsDummy() || block.uid != 4);
  }
  EXPECT_EQ(realCount, 8UL);
}

TEST(CircuitORAM, ForcedLegacyAndOneTimeHashBatchStashModes) {
  constexpr uint64_t elementCount = 128;
#ifdef DISK_IO
  if (EM::Backend::g_DefaultBackend != nullptr) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(1UL << 24);
  ODSL::CircuitORAM::ORAM<uint64_t> oram(elementCount, 1UL << 10);
#else
  ODSL::CircuitORAM::ORAM<uint64_t> oram(elementCount);
#endif
  std::vector<uint64_t> positions(elementCount);
  for (uint64_t uid = 0; uid < elementCount; ++uid) {
    positions[uid] = oram.Write(uid, 10'000 + uid);
  }

  auto runBatch = [&](uint64_t begin,
                      ODSL::CircuitORAM::BatchStashAccessMode mode) {
    constexpr uint64_t batchSize = 64;
    std::array<uint64_t, batchSize> uids;
    std::array<uint64_t, batchSize> batchPositions;
    std::array<uint64_t, batchSize> newPositions;
    std::array<uint64_t, batchSize> values{};
    for (uint64_t i = 0; i < batchSize; ++i) {
      uids[i] = begin + i;
      batchPositions[i] = positions[uids[i]];
    }

    oram.BatchReadAndRemove(batchSize, batchPositions.data(), uids.data(),
                            values.data(), mode);
    for (uint64_t i = 0; i < batchSize; ++i) {
      EXPECT_EQ(values[i], 10'000 + uids[i]);
    }
    oram.GetRandNewPoses(newPositions.data(), batchSize);
    oram.BatchWriteBack(batchSize, uids.data(), newPositions.data(),
                        values.data(), std::vector<bool>(batchSize, true));
    for (uint64_t i = 0; i < batchSize; ++i) {
      positions[uids[i]] = newPositions[i];
    }
  };

  runBatch(0, ODSL::CircuitORAM::BatchStashAccessMode::LegacyScan);
  runBatch(64, ODSL::CircuitORAM::BatchStashAccessMode::OneTimeHash);

  for (uint64_t uid = 0; uid < elementCount; ++uid) {
    uint64_t value = 0;
    positions[uid] = oram.Read(positions[uid], uid, value);
    EXPECT_EQ(value, 10'000 + uid);
  }

  // A forced empty hash batch still performs a well-formed build/restore and
  // must not require non-null request buffers.
  oram.BatchReadAndRemove(0, nullptr, nullptr, nullptr,
                          ODSL::CircuitORAM::BatchStashAccessMode::OneTimeHash);
}

}  // namespace
