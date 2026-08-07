#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "odsl/circuit_oram.hpp"
#include "odsl/one_time_stash_hash.hpp"

namespace {

using Clock = std::chrono::steady_clock;

uint64_t benchmarkSink = 0;

struct Options {
  size_t rounds = 100;
  size_t stashLoad = 23;
  size_t payloadBytes = 0;
  size_t batchSize = 0;
  bool micro = true;
  bool endToEnd = false;
  size_t oramRounds = 10;
};

template <size_t Words>
struct Payload {
  std::array<uint64_t, Words> words{};

  static consteval Payload DUMMY() { return {}; }
};

template <typename Setup, typename Func>
double MeasureMicroseconds(size_t rounds, Setup&& setup, Func&& func) {
  std::chrono::nanoseconds elapsed{};
  for (size_t round = 0; round < rounds; ++round) {
    setup();
    const auto start = Clock::now();
    const uint64_t result = func();
    const auto end = Clock::now();
    elapsed += end - start;
    benchmarkSink += result;
  }
  return std::chrono::duration<double, std::micro>(elapsed).count() / rounds;
}

size_t RoundsForBatch(size_t baseRounds, size_t batchSize) {
  return std::max<size_t>(3,
                          baseRounds * 256 / std::max<size_t>(256, batchSize));
}

template <size_t Words>
void RunMicrobenchmarkForPayload(const Options& options) {
  constexpr size_t stashCapacity = 46;
  using Block = ODSL::Block<Payload<Words>, uint64_t, uint64_t>;
  using Hash =
      ODSL::OneTimeStashHash<Block, uint64_t, stashCapacity, 4, 80, 4, 1>;

  const size_t stashLoad = std::min(options.stashLoad, stashCapacity);
  std::array<Block, stashCapacity> initialStash{};
  for (size_t i = 0; i < stashLoad; ++i) {
    Payload<Words> payload;
    payload.words[0] = i + 1;
    initialStash[i] = Block(payload, i, i);
  }

  for (bool includeHits : {false, true}) {
    for (size_t batchSize :
         {1UL, 4UL, 16UL, 64UL, 256UL, 1024UL, 4096UL, 16384UL}) {
      if (options.batchSize != 0 && options.batchSize != batchSize) {
        continue;
      }
      std::vector<uint64_t> queries(batchSize);
      const size_t hitCount = includeHits ? std::min(stashLoad, batchSize) : 0;
      for (size_t i = 0; i < batchSize; ++i) {
        queries[i] = i < hitCount ? i : UINT64_C(1000000) + i;
      }

      const size_t rounds = RoundsForBatch(options.rounds, batchSize);
      std::array<Block, stashCapacity> legacyStash{};
      const double legacyUs = MeasureMicroseconds(
          rounds, [&]() { legacyStash = initialStash; },
          [&]() {
            uint64_t foundCount = 0;
            Payload<Words> out;
            for (uint64_t uid : queries) {
              foundCount += ODSL::ReadElementAndRemoveFromPath(
                  legacyStash.begin(), legacyStash.end(), uid, out);
            }
            return foundCount + out.words[0];
          });

      std::array<Block, stashCapacity> hashStash{};
      size_t fallbackCount = 0;
      const double hashUs = MeasureMicroseconds(
          rounds, [&]() { hashStash = initialStash; },
          [&]() {
            uint64_t foundCount = 0;
            Payload<Words> out;
            Hash hash;
            if (hash.Build(hashStash.data())) {
              hash.PrepareBatch(queries.data(), queries.size());
              for (size_t i = 0; i < queries.size(); ++i) {
                foundCount +=
                    hash.ReadAndRemovePrepared(queries[i], i, true, out);
              }
              hash.Restore();
            } else {
              ++fallbackCount;
              for (uint64_t uid : queries) {
                foundCount += ODSL::ReadElementAndRemoveFromPath(
                    hashStash.begin(), hashStash.end(), uid, out);
              }
            }
            return foundCount + out.words[0];
          });

      std::cout << "micro," << sizeof(Payload<Words>) << ',' << stashCapacity
                << ',' << stashLoad << ','
                << (includeHits ? "hits_then_misses" : "misses") << ','
                << batchSize << ',' << rounds << ',' << std::fixed
                << std::setprecision(3) << legacyUs << ',' << hashUs << ','
                << hashUs / legacyUs << ',' << fallbackCount << '\n';
    }
  }
}

struct OramTimings {
  double legacyUs;
  double hashUs;
};

template <size_t Words>
OramTimings RunOramReadBenchmark(size_t batchSize, size_t rounds) {
  constexpr size_t stashCapacity = 46;
  constexpr size_t elementCount = 32768;
  using Oram = ODSL::CircuitORAM::ORAM<Payload<Words>, 2, stashCapacity,
                                       uint64_t, uint64_t, 4096, false>;

#ifdef DISK_IO
  if (EM::Backend::g_DefaultBackend != nullptr) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(1UL << 28);
  Oram oram(elementCount, 1UL << 20);
#else
  Oram oram(elementCount);
#endif
  std::vector<uint64_t> positions(elementCount);
  for (uint64_t uid = 0; uid < elementCount; ++uid) {
    Payload<Words> payload;
    payload.words[0] = uid + 1;
    positions[uid] = oram.Write(uid, payload);
  }

  std::vector<uint64_t> uids(batchSize);
  std::vector<uint64_t> batchPositions(batchSize);
  std::vector<uint64_t> newPositions(batchSize);
  std::vector<Payload<Words>> values(batchSize);
  const std::vector<bool> writeBackFlags(batchSize, true);
  std::chrono::nanoseconds legacyElapsed{};
  std::chrono::nanoseconds hashElapsed{};

  for (size_t round = 0; round < rounds; ++round) {
    const size_t begin = (round * batchSize) % (elementCount - batchSize + 1);
    for (size_t i = 0; i < batchSize; ++i) {
      uids[i] = begin + i;
    }
    const bool hashFirst = round & 1;
    for (size_t pass = 0; pass < 2; ++pass) {
      const bool useHash = (pass == 0) == hashFirst;
      const auto mode =
          useHash ? ODSL::CircuitORAM::BatchStashAccessMode::OneTimeHash
                  : ODSL::CircuitORAM::BatchStashAccessMode::LegacyScan;
      for (size_t i = 0; i < batchSize; ++i) {
        batchPositions[i] = positions[uids[i]];
      }
      oram.GetRandNewPoses(newPositions.data(), batchSize);

      const auto start = Clock::now();
      oram.BatchReadAndRemove(batchSize, batchPositions.data(), uids.data(),
                              values.data(), mode);
      const auto end = Clock::now();
      if (useHash) {
        hashElapsed += end - start;
      } else {
        legacyElapsed += end - start;
      }

      oram.BatchWriteBack(batchSize, uids.data(), newPositions.data(),
                          values.data(), writeBackFlags);
      for (size_t i = 0; i < batchSize; ++i) {
        positions[uids[i]] = newPositions[i];
        benchmarkSink += values[i].words[0];
      }
    }
  }

  return {
      std::chrono::duration<double, std::micro>(legacyElapsed).count() / rounds,
      std::chrono::duration<double, std::micro>(hashElapsed).count() / rounds};
}

template <size_t Words>
void RunEndToEndForPayload(const Options& options) {
  for (size_t batchSize :
       {64UL, 256UL, 1024UL, 1536UL, 2048UL, 3072UL, 4096UL, 16384UL}) {
    if (options.batchSize != 0 && options.batchSize != batchSize) {
      continue;
    }
    const size_t rounds =
        std::max<size_t>(10, options.oramRounds * 64 / batchSize);
    const OramTimings timings = RunOramReadBenchmark<Words>(batchSize, rounds);
    std::cout << "oram_read," << sizeof(Payload<Words>)
              << ",46,-,existing_n32768," << batchSize << ',' << rounds << ','
              << std::fixed << std::setprecision(3) << timings.legacyUs << ','
              << timings.hashUs << ',' << timings.hashUs / timings.legacyUs
              << ",0\n";
  }
}

size_t ParseSize(std::string_view value, std::string_view option) {
  std::string text(value);
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || parsed == 0) {
    throw std::invalid_argument("invalid value for " + std::string(option));
  }
  return static_cast<size_t>(parsed);
}

Options ParseOptions(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if (arg == "--end-to-end") {
      options.endToEnd = true;
    } else if (arg == "--end-to-end-only") {
      options.micro = false;
      options.endToEnd = true;
    } else if (arg.starts_with("--rounds=")) {
      options.rounds = ParseSize(arg.substr(9), "--rounds");
    } else if (arg.starts_with("--oram-rounds=")) {
      options.oramRounds = ParseSize(arg.substr(14), "--oram-rounds");
    } else if (arg.starts_with("--stash-load=")) {
      options.stashLoad = ParseSize(arg.substr(13), "--stash-load");
    } else if (arg.starts_with("--payload-bytes=")) {
      options.payloadBytes = ParseSize(arg.substr(16), "--payload-bytes");
    } else if (arg.starts_with("--batch-size=")) {
      options.batchSize = ParseSize(arg.substr(13), "--batch-size");
    } else {
      throw std::invalid_argument("unknown option: " + std::string(arg));
    }
  }
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = ParseOptions(argc, argv);
    std::cout
        << "scope,payload_bytes,stash_capacity,stash_load,pattern,batch_size,"
           "rounds,legacy_us,hash_us,hash_over_legacy,hash_fallbacks\n";
    if (options.micro) {
      if (options.payloadBytes == 0 || options.payloadBytes == 8) {
        RunMicrobenchmarkForPayload<1>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 64) {
        RunMicrobenchmarkForPayload<8>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 256) {
        RunMicrobenchmarkForPayload<32>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 1024) {
        RunMicrobenchmarkForPayload<128>(options);
      }
    }
    if (options.endToEnd) {
      if (options.payloadBytes == 0 || options.payloadBytes == 8) {
        RunEndToEndForPayload<1>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 64) {
        RunEndToEndForPayload<8>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 128) {
        RunEndToEndForPayload<16>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 256) {
        RunEndToEndForPayload<32>(options);
      }
      if (options.payloadBytes == 0 || options.payloadBytes == 1024) {
        RunEndToEndForPayload<128>(options);
      }
    }
    std::cerr << "benchmark_sink=" << benchmarkSink << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
