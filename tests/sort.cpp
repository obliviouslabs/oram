#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/ca_bucket_sort.hpp"
#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/algorithm/kway_distri_sort.hpp"
#include "testutils.hpp"

using namespace EM::Algorithm;
using namespace EM::NonCachedVector;
using namespace std;

TEST(TestSort, Compact) {
  uint64_t N = 7;
  uint64_t num1 = 5;
  vector<uint64_t> vec(N, 0);
  for (size_t i = 0; i < num1; ++i) {
    vec[i] = 1;
  }
  fisherYatesShuffle(vec.begin(), vec.end());
  const auto isMarked = [=](const auto& element) { return element; };
  for (auto n : vec) {
    cout << n << " ";
  }
  cout << endl;
  GoodrichCompact(vec.begin(), vec.end(), isMarked);
  for (auto n : vec) {
    cout << n << " ";
  }
  cout << endl;
}

TEST(TestSort, MergeSplit) {
  for (uint64_t Z = 400; Z < 800; ++Z) {
    const uint64_t N = 2 * Z;
    const uint64_t beginIdxLeft = 0;
    const uint64_t beginIdxRight = Z;
    for (uint64_t bitMask = 1; bitMask <= 8; bitMask *= 2) {
      TaggedT<uint64_t> defaultval;
      vector<TaggedT<uint64_t>> v(N, defaultval);
      unordered_map<uint64_t, long long> tagCount;
      long long dummyCount = 0;
      for (uint64_t i = 0; i < 2 * Z; ++i) {
        uint64_t idx = i < Z ? i + beginIdxLeft : i - Z + beginIdxRight;
        uint64_t rnd = UniformRandom();

        v[idx].setTag(rnd);
        uint64_t randVal = v[idx].tag * 3;
        if (randVal % 929 <= 929 / 2) {
          v[idx].v = randVal;
          ++tagCount[v[idx].tag];
        } else {
          v[idx].setDummy();
          ++dummyCount;
        }
      }

      fisherYatesShuffle(v.begin(), v.end());
      MergeSplitInPlace(v.begin(), v.end(), bitMask);

      for (uint64_t i = 0; i < 2 * Z; ++i) {
        uint64_t idx = i < Z ? i + beginIdxLeft : i - Z + beginIdxRight;
        if (v[idx].isDummy()) {
          ASSERT_GE(--dummyCount, 0);
        } else {
          ASSERT_GE(--tagCount[v[idx].tag], 0);
        }
      }

      // check partitioned and tags/values unchanged
      for (uint64_t i = 0; i < Z; ++i) {
        uint64_t idx = i + beginIdxLeft;
        ASSERT_TRUE((v[idx].tag & bitMask) == 0);
        if (!v[idx].isDummy()) {
          ASSERT_TRUE(v[idx].tag * 3 == v[idx].v);
        }
      }
      for (uint64_t i = 0; i < Z; ++i) {
        uint64_t idx = i + beginIdxRight;
        ASSERT_TRUE((v[idx].tag & bitMask) != 0);
        if (!v[idx].isDummy()) {
          ASSERT_TRUE(v[idx].tag * 3 == v[idx].v);
        }
      }
    }
  }
}

TEST(TestSort, MergeSplitUsingPivots) {
  for (uint64_t Z = 400; Z < 800; ++Z) {
    const uint64_t N = 2 * Z;
    const uint64_t beginIdxLeft = 0;
    const uint64_t beginIdxRight = Z;
    for (uint64_t bitMask = 1; bitMask <= 8; bitMask *= 2) {
      Block<uint64_t> defaultval;
      vector<Block<uint64_t>> v(N, defaultval);
      unordered_map<uint64_t, long long> tagCount;
      long long dummyCount = 0;
      for (uint64_t i = 0; i < 2 * Z; ++i) {
        uint64_t idx = i < Z ? i + beginIdxLeft : i - Z + beginIdxRight;
        uint64_t rnd = UniformRandom(200);

        if (UniformRandom() % 929 <= 929 / 2) {
          v[idx].setData(rnd);
          ++tagCount[rnd];
        } else {
          v[idx].setDummy();
          ++dummyCount;
        }
      }

      Block<uint64_t> pivot;
      pivot.setData(100UL);

      fisherYatesShuffle(v.begin(), v.end());

      MergeSplitInPlace(v.begin(), v.end(), pivot);
      // check it's a permutation
      for (uint64_t i = 0; i < 2 * Z; ++i) {
        uint64_t idx = i < Z ? i + beginIdxLeft : i - Z + beginIdxRight;
        if (v[idx].isDummy()) {
          ASSERT_GE(--dummyCount, 0);
        } else {
          ASSERT_GE(--tagCount[v[idx].data], 0);
        }
      }

      // check partitioned
      for (uint64_t i = 0; i < Z; ++i) {
        uint64_t idx = i + beginIdxLeft;
        if (!v[idx].isDummy()) ASSERT_TRUE(v[idx] < pivot);
      }
      for (uint64_t i = 0; i < Z; ++i) {
        uint64_t idx = i + beginIdxRight;
        if (!v[idx].isDummy()) ASSERT_FALSE(v[idx] < pivot);
      }
    }
  }
}

TEST(TestSort, TestButterflyWaySolver) {
  ButterflyWaySolver solver1({3, 4, 5, 6, 7, 8}, {8, 9, 10, 11, 12, 13}, 0,
                             100000);
  vector<vector<uint64_t>> result;
  solver1.solve(result, 510);
  for (auto sub : result) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << endl;
  ButterflyWaySolver solver2({3, 4, 5, 6, 7, 8}, {8, 9, 10, 11, 12, 13}, 40,
                             30);
  solver2.solve(result, 510);
  for (auto sub : result) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << endl;
  solver2.solve(result, 2001);
  for (auto sub : result) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << endl;
  solver2.solve(result, 3);
  for (auto sub : result) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << endl;
  ButterflyWaySolver solver3({3, 4, 5, 6, 7, 8}, {8, 9, 10, 11, 12, 13}, 1, 30);
  solver3.solve(result, 510);
  for (auto sub : result) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << endl;
}

TEST(TestSort, TestKWayButterflyParams) {
  KWayButterflyParams params = bestKWayButterflyParams(800UL << 19);
  for (auto sub : params.ways) {
    for (auto way : sub) {
      cout << way << " ";
    }
    cout << endl;
  }
  cout << "Z=" << params.Z << endl;
}

TEST(TestSort, EulerPath) {
  // single loop
  EdgeRec rec(4);
  rec.flipEdge(0, 1);
  rec.flipEdge(2, 3);
  rec.flipEdge(2, 1);
  rec.flipEdge(0, 3);
  rec.print();
  EdgeRec path = rec.EulerPath(4);
  rec.printPath(path);
  // two loops
  EdgeRec rec2(6);
  rec2.flipEdge(0, 1);
  rec2.flipEdge(1, 3);
  rec2.flipEdge(0, 3);
  rec2.flipEdge(2, 4);
  rec2.flipEdge(2, 5);
  rec2.flipEdge(4, 5);
  rec2.print();
  path = rec2.EulerPath(6);
  rec2.printPath(path);
  // star
  EdgeRec rec3(5);
  rec3.flipEdge(0, 2);
  rec3.flipEdge(2, 4);
  rec3.flipEdge(4, 1);
  rec3.flipEdge(1, 3);
  rec3.flipEdge(3, 0);
  rec3.flipEdge(0, 1);
  rec3.flipEdge(1, 2);
  rec3.flipEdge(2, 3);
  rec3.flipEdge(3, 4);
  rec3.flipEdge(4, 0);
  rec3.print();
  path = rec3.EulerPath(10);
  rec3.printPath(path);
}

TEST(TestSort, KWayInterleaveSepMarks) {
  vector<uint8_t> marks = {3, 3, 1, 4, 3, 0, 2, 0, 4, 1,
                           0, 3, 1, 2, 4, 4, 0, 1, 2, 2};
  vector<uint64_t> v(marks.size());
  for (size_t i = 0; i < marks.size(); ++i) {
    v[i] = (uint64_t)marks[i] * 10UL;
  }

  Interleave(v.begin(), v.end(), marks.begin(), marks.end(), 5);
  for (auto ele : v) {
    cout << ele << " ";
  }
  cout << endl;
  for (size_t way = 3; way <= 8; ++way) {
    vector<uint64_t> vLarge(16384 * way);
    vector<uint64_t> markLarge(16384 * way);
    for (int round = 0; round < 50; ++round) {
      for (size_t i = 0; i < way; ++i) {
        for (size_t j = 0; j < 16384; ++j) {
          markLarge[i * 16384 + j] = i;
        }
      }
      fisherYatesShuffle(markLarge.begin(), markLarge.end());
      for (size_t i = 0; i < markLarge.size(); ++i) {
        vLarge[i] = markLarge[i] * 10;
      }
      Interleave(vLarge.begin(), vLarge.end(), markLarge.begin(),
                 markLarge.end(), way);
      for (uint i = 0; i < vLarge.size(); ++i) {
        ASSERT_EQ(vLarge[i] / 10, i % way);
      }
    }
  }
}

TEST(TestSort, KWayMergeSplit) {
  for (size_t way = 3; way <= 8; ++way) {
    vector<TaggedT<uint64_t>> vLarge(16384 * way);
    for (int round = 0; round < 50; ++round) {
      for (size_t i = 0; i < way; ++i) {
        for (size_t j = 0; j < 16384; ++j) {
          vLarge[i * 16384 + j].setData(i);
          vLarge[i * 16384 + j].setTag(i + way * 1000);
          if (UniformRandom() % 4 == 0) {
            vLarge[i * 16384 + j].setDummy();
          }
        }
      }
      fisherYatesShuffle(vLarge.begin(), vLarge.end());
      vector<TaggedT<uint64_t>>::iterator begins[8];

      for (size_t i = 0; i < way; ++i) {
        begins[i] = vLarge.begin() + i * 16384;
      }
      MergeSplitKWay(begins, way, 16384);
      for (size_t i = 0; i < vLarge.size(); ++i) {
        if (!vLarge[i].isDummy()) {
          ASSERT_EQ(vLarge[i].v, i / 16384);
        }
      }
    }
    printf("way %zu passes\n", way);
  }
}

void testKWayParams() {
  cout.precision(4);
  size_t b = 128;
  for (int lgSize = 6; lgSize <= 10; ++lgSize) {
    size_t N = pow(10, lgSize);
    double target = -60;
    size_t M = 0x7500000 / (b + 8);

    KWayButterflyParams params;
    params = bestKWayButterflyParams(N, M, b, target);

    size_t Z = params.Z;
    double padRate = double(params.totalBucket * params.Z) / N - 1.0;
    auto satisfy = [=](size_t bucketCount) {
      return failureProbButterflySort(Z, N, bucketCount) < target;
    };

    size_t minBucketCount = lowerBound(N / Z + 1, 3 * N / Z, satisfy);
    double minPadRate = (double)minBucketCount * Z / N - 1;
    size_t Zr = ceil(Z / (1 + padRate));
    double failProb = failureProbButterflySort(params.Z, N, params.totalBucket);
    cout << "$10^{" << lgSize << "}$ & " << params.Z << " & " << minPadRate
         << " & " << padRate << " & ";
    bool init = true;
    for (auto sub : params.ways) {
      if (!init) {
        cout << ")\\times";
      } else {
        cout << "$";
      }
      init = false;
      cout << "(";
      bool inInit = true;
      for (auto way : sub) {
        if (!inInit) {
          cout << "\\times";
        }
        inInit = false;
        cout << way;
      }
    }
    cout << ")$";
    cout << endl;
  }
}

void testDistriParams() {
  cout.precision(4);
  size_t b = 128;
  for (size_t b = 128; b <= 128; b = b * 3 / 2) {
    size_t B = 4096 / b;
    for (int lgSize = 6; lgSize <= 10; ++lgSize) {
      size_t N = pow(10, lgSize);
      double target = -60;
      size_t wrappedSize = b % 32 == 16 ? b + 16 : b + 8;
      size_t M = 0x7500000 / wrappedSize;

      DistriParams params = bestDistriParams(N, M, B, b, target);
      size_t Z = params.Z;
      double padRate = double(params.totalBucket * params.Z) / N - 1.0;

      cout << "$10^{" << lgSize << "}$ & " << params.Z << " & "
           << params.samplingRatio << " & " << padRate << " & ";
      bool init = true;
      for (auto sub : params.ways) {
        if (!init) {
          cout << ")\\times";
        } else {
          cout << "$";
        }
        init = false;
        cout << "(";
        bool inInit = true;
        for (auto way : sub) {
          if (!inInit) {
            cout << "\\times";
          }
          inInit = false;
          cout << way;
        }
      }
      cout << ")\\times";
      cout << params.totalBucket / params.totalPartition;
      cout << "$" << endl;
    }
  }
}

TEST(TestSort, KWayParams) {
  printf("--------Butterfly-----------\n");
  testKWayParams();
  printf("--------Distri-----------\n");
  testDistriParams();
}

TEST(TestSort, testAuth) {
  struct Page {
    uint64_t pages[128];
    using Encrypted_t = FreshEncrypted<Page>;
  };
  EM::MemoryServer::NonCachedServerFrontendInstance<
      Page, ::EM::Backend::MemServerBackend, true, true>
      front(*EM::Backend::g_DefaultBackend, 4);
  Page page;
  for (size_t i = 0; i < 128; ++i) {
    page.pages[i] = i;
  }
  front.Write(1, page, 0);
  front.flushWrite();
  for (size_t i = 0; i < 12; ++i) {
    cout << (int)front.nounce.bytes[i] << " ";
  }
  cout << endl;
  front.Read(1, page, 0);
  front.flushRead();
  for (size_t i = 0; i < 12; ++i) {
    cout << (int)front.nounce.bytes[i] << " ";
  }
  cout << endl;
  for (size_t i = 0; i < 128; ++i) {
    cout << page.pages[i] << " ";
  }
  cout << endl;
}

TEST(TestSort, MergePermute) {
  std::vector<int64_t> data;
  std::vector<uint8_t> marks;
  for (int i = 2; i <= 8; ++i) {
    for (int round = 0; round < 1000; ++round) {
      for (int j = 0; j < i; ++j) {
        marks.push_back(j);
      }
      fisherYatesShuffle(marks.begin(), marks.end());
      for (auto mark : marks) {
        data.push_back(mark * 10);
      }
      Permute(data.begin(), data.end(), marks.begin(), marks.end());
      for (int j = 0; j < i; ++j) {
        ASSERT_EQ(data[j], j * 10);
      }
      marks.clear();
      data.clear();
    }
  }
}