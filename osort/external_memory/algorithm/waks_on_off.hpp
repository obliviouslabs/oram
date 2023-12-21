#pragma once
#include <unordered_map>
#include <vector>

#include "bitonic.hpp"
#include "external_memory/stdvector.hpp"
#include "sort_def.hpp"

namespace EM::Algorithm {

class WaksOnOff {
 public:
  struct ControlCircuit {
    std::vector<bool> in;
    std::vector<bool> out;
    const ControlCircuit* top = NULL;
    const ControlCircuit* bot = NULL;
    ~ControlCircuit() {
      if (top != NULL) {
        delete top;
      }
      if (bot != NULL) {
        delete bot;
      }
    }
  };

  static const ControlCircuit* WaksShuffleOffline(uint64_t n) {
    StdVector<uint64_t> P(n);
    for (uint64_t i = 0; i < P.size(); ++i) {
      P[i] = i;
    }
    OrShuffle(P.begin(), P.end());
    return WaksShuffleOffline(P.begin(), P.end());
  }

  template <typename Iterator>
  static const ControlCircuit* WaksShuffleOffline(Iterator permBegin,
                                                  Iterator permEnd) {
    const ControlCircuit* C = ControlBits(permBegin, permEnd);
    return C;
  }

  template <typename Iterator>
  static void WaksShuffleOnline(Iterator begin, Iterator end,
                                const ControlCircuit* C) {
    uint64_t n = end - begin;
    uint64_t k = divRoundUp(n, 2);
    for (uint64_t i = 0; i < k - 1; ++i) {
      condSwap(C->in[i], *(begin + i), *(begin + (k + i)));
    }
    if (n > 2) {
      WaksShuffleOnline(begin, begin + k, C->top);
      WaksShuffleOnline(begin + k, end, C->bot);
    }
    for (uint64_t i = 0; i < n - k; ++i) {
      condSwap(C->out[i], *(begin + i), *(begin + (k + i)));
    }
  }

  template <typename Iterator>
  static void WaksOnOffShuffle(Iterator begin, Iterator end) {
    const ControlCircuit* C = WaksShuffleOffline(end - begin);
    WaksShuffleOnline(begin, end, C);
    delete C;
  }

  static void printCircuit(const ControlCircuit* C, int depth = 0) {
    if (!C) {
      printf("{}");
      return;
    }
    printf("{\"in\":[");
    for (size_t i = 0; i < C->in.size(); ++i) {
      if (i == C->in.size() - 1) {
        printf("%d", C->in[i]);
      } else {
        printf("%d, ", C->in[i]);
      }
    }
    printf("],");
    printf("\n");
    printf("\"top\":");
    printCircuit(C->top, depth + 1);
    printf(",\n");
    printf("\"bot\":");
    printCircuit(C->bot, depth + 1);
    printf(",\n\"out\": [");
    for (size_t i = 0; i < C->out.size(); ++i) {
      if (i == C->out.size() - 1) {
        printf("%d", C->out[i]);
      } else {
        printf("%d, ", C->out[i]);
      }
    }
    printf("]}\n");
  }

 private:
  /**
   * Compute Waksman control bits for permutation P at recursion level d. i maps
   * to a value encoded in P[i]. Initially call with d = 0. Modifies P.
   * Oblivious to values in P but not to |P | or d.
   */
  template <typename Iterator>
  static const ControlCircuit* ControlBits(Iterator PBegin, Iterator PEnd,
                                           uint64_t d = 0) {
    uint64_t n = PEnd - PBegin;
    uint64_t k = divRoundUp(n, 2);
    if (n < 2) {
      return NULL;
    }
    ControlCircuit* circuit = new ControlCircuit();
    // Compute control bits for input-layer switches
    circuit->in = InBits(PBegin, PEnd, d);
    // Apply input-layer switches to P
    if (n > 2) {
      for (uint64_t i = 0; i < k - 1; i++) {
        condSwap(circuit->in[i], *(PBegin + i), *(PBegin + (k + i)));
      }
    }
    // Reduce P values modulo k, storing bits to undo later
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t b = (*(PBegin + i) >> d) >= k;
      *(PBegin + i) = ((*(PBegin + i) - (OSelect(b, 0, k) << d)) << 1) | b;
    }
    // Compute control bits for top and bottom subnetworks
    if (n > 2) {
      circuit->top = ControlBits(PBegin, PBegin + k, d + 1);
      circuit->bot = ControlBits(PBegin + k, PEnd, d + 1);
    }
    circuit->out.resize(n - k);
    // Compute control bits for output-layer switches. Create array Cout of
    // length n − k. Apply output-layer switches to P.
    for (uint64_t i = 0; i < n - k; ++i) {
      circuit->out[i] = *(PBegin + i) & 1;
      condSwap(circuit->out[i], *(PBegin + i), *(PBegin + (k + i)));
    }
    // Undo reduction of P values modulo k
    for (uint64_t i = 0; i < n; ++i) {
      bool b = (*(PBegin + i) & 1);
      *(PBegin + i) = (*(PBegin + i) >> 1) + (OSelect(b, 0, k) << d);
    }
    return circuit;
  }

  typedef std::vector<TaggedT<std::pair<uint64_t, uint64_t>>>
      ForwardLookUpStruct;

  typedef std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>>
      ReverseLookUpStruct;

  typedef std::vector<uint64_t> UnselectedCountsStruct;

  template <typename Iterator>
  static std::vector<bool> InBits(Iterator PBegin, Iterator PEnd, uint64_t d) {
    uint64_t n = PEnd - PBegin;
    uint64_t k = divRoundUp(n, 2);
    std::vector<bool> C(k - 1);
    if (n <= 2) {
      return C;
    }
    // Switch numbers of control bits
    PRP<uint64_t> forwardPrp, backwardPrp;
    std::vector<uint64_t> S(k - 1);
    ForwardLookUpStruct F = CreateForwardLookUp(PBegin, PEnd, d, forwardPrp);
    ReverseLookUpStruct R = CreateReverseLookUp(F, backwardPrp);
    // printf("forward lookup\n");
    // for (auto td : F) {
    //   printf("%lu (%lu, %lu)\n", td.tag, td.v.first, td.v.second);
    // }
    // printf("reverse lookup\n");
    // for (auto [k, v] : R) {
    //   printf("%lu (%lu, %lu)\n", k, v.first, v.second);
    // }
    UnselectedCountsStruct U = CreateUnselectedCounts(2 * k);
    uint64_t l = ForwardOrRand(F, U, k - 1, 0, forwardPrp);  // map k-1
    uint64_t f = F[l].v.first;
    uint64_t g = F[l].v.second;

    DecUnselectedCounts(U.begin(), U.end(), l);  // mark k-1 -> g as selected
    uint64_t c = f;  // cycle start is f (which contains k-1)
    uint64_t r = g + OSelect(g < k, -k, k);       // r is associate of g
    auto [s, l2] = R[backwardPrp.prp(r)];         // reverse map of r
    DecUnselectedCounts(U.begin(), U.end(), l2);  // s->r map as selected
    f = s + OSelect(s < k, -k, k);                // f is associate of s
    // printf("f updated = %lu\n", f);
    // repeat back and forth process to compute the control bits
    for (uint64_t i = 0; i < k - 1; ++i) {
      bool b = f == c;
      uint64_t l = ForwardOrRand(F, U, f, b, forwardPrp);
      f = F[l].v.first;
      uint64_t g = F[l].v.second;
      DecUnselectedCounts(U.begin(), U.end(), l);
      c = OSelect(b, c, f);
      uint64_t r = g + OSelect(g < k, -k, k);
      auto [s, l2] = R[backwardPrp.prp(r)];
      DecUnselectedCounts(U.begin(), U.end(), l2);
      C[i] = f >= k;
      S[i] = f - OSelect(f >= k, 0, k);
      f = s + OSelect(s < k, -k, k);
      // printf("complete dec unselected counts %ld time\n", i);
    }
    std::vector<TaggedT<bool>> CTagged(k - 1);
    for (uint64_t i = 0; i < k - 1; ++i) {
      CTagged[i].tag = S[i];
      CTagged[i].v = C[i];
    }
    BitonicSort(CTagged.begin(), CTagged.end(),
                [](const auto& a, const auto& b) { return a.tag < b.tag; });
    for (uint64_t i = 0; i < k - 1; ++i) {
      C[i] = CTagged[i].v;
    }
    return C;
  }

  template <typename Iterator>
  static ForwardLookUpStruct CreateForwardLookUp(
      Iterator PBegin, Iterator PEnd, uint64_t d,
      const PRP<uint64_t>& forwardPrp) {
    uint64_t n = PEnd - PBegin;
    uint64_t k = divRoundUp(n, 2);
    uint64_t m = 2 * k;
    ForwardLookUpStruct forwardLookUp(m);

    for (uint64_t i = 0; i < m; ++i) {
      forwardLookUp[i].tag = forwardPrp.prp(i);
    }
    for (uint64_t i = 0; i < n; ++i) {
      forwardLookUp[i].v = std::make_pair(i, *(PBegin + i) >> d);
    }
    if (n != m) {
      forwardLookUp[n].v = std::make_pair(n, n);
    }
    auto cmpTag = [](const auto& a, const auto& b) { return a.tag < b.tag; };
    BitonicSort(forwardLookUp.begin(), forwardLookUp.end(), cmpTag);
    return forwardLookUp;
  }

  static ReverseLookUpStruct CreateReverseLookUp(
      ForwardLookUpStruct& forwardLookUp, const PRP<uint64_t>& reversePrp) {
    ReverseLookUpStruct reverseLookUp;
    uint64_t n = forwardLookUp.size();
    reverseLookUp.reserve(n);
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t x = forwardLookUp[i].v.first;
      uint64_t y = forwardLookUp[i].v.second;
      uint64_t z = reversePrp.prp(y);

      reverseLookUp[z] = std::make_pair(x, i);
    }
    return reverseLookUp;
  }

  static UnselectedCountsStruct CreateUnselectedCounts(uint64_t n) {
    UnselectedCountsStruct U(n);
    U[n - 1] = n;
    if (n < 2) {
      return U;
    }
    CreateUnselectedCountsHelper(U.begin(), U.end());
    return U;
  }

  template <typename Iterator>
  static void CreateUnselectedCountsHelper(Iterator begin, Iterator end) {
    uint64_t n = end - begin;
    uint64_t k = divRoundUp(n, 2);
    if (n < 2) {
      return;
    }
    *(begin + k - 1) = k;
    CreateUnselectedCountsHelper(begin, begin + k);
    CreateUnselectedCountsHelper(begin + k, end);
  }

  /**
   *  If b = 0, look up for-
 ward map from f under permutation encoded in F. Else (i.e., if b = 1), look up
 a uniformly random mapping in F from among those that U indicates are
 unselected. Returns mapping index l.
  */
  static uint64_t ForwardOrRand(const ForwardLookUpStruct& F,
                                const UnselectedCountsStruct& U, uint64_t f,
                                bool b, const PRP<uint64_t>& forwardPrp) {
    // printf("forward or rand with b = %d and f = %lu\n", b, f);
    uint64_t n = F.size();
    uint64_t h = forwardPrp.prp(f);
    // printf("h = %lu\n", h);
    uint64_t i = 0;
    uint64_t j = n - 1;
    uint64_t u = U[n - 1];
    uint64_t rho = UniformRandom(u - 1);
    uint64_t l;
    while (true) {
      l = (i + j) / 2;
      if (i == j) {
        break;
      }
      bool oSelectFalse = F[l].tag < h;
      bool oSelectTrue = U[l] <= rho;
      bool z = OSelect(b, oSelectFalse, oSelectTrue);
      if (!z) {
        j = l;
        u = U[l];
      } else {
        i = l + 1;
        u = u - U[l];
        rho = rho - U[l];
      }
    }
    return l;
  }

  template <typename Iterator>
  static void DecUnselectedCounts(Iterator begin, Iterator end, uint64_t l) {
    Assert(begin != end);
    uint64_t n = end - begin;
    Assert(l < n);
    uint64_t k = divRoundUp(n, 2);
    if (l < k) {
      --*(end - 1);
      if (n > 1) {
        DecUnselectedCounts(begin, begin + k, l);
      }
    } else {
      DecUnselectedCounts(begin + k, end, l - k);
    }
  }

  static uint64_t OSelect(bool b, uint64_t x, uint64_t y) { return b ? y : x; }
};

template <typename Iterator>
void WaksOnOffShuffle(Iterator begin, Iterator end) {
  WaksOnOff::WaksOnOffShuffle(begin, end);
}

template <typename Vec>
void WaksOnOffShuffle(Vec& v) {
  WaksOnOff::WaksOnOffShuffle(v.begin(), v.end());
}
}  // namespace EM::Algorithm