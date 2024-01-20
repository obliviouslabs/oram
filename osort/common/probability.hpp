#pragma once

#include "utils.hpp"

double UniformRandomDouble() {
  return (double)UniformRandom() / (double)UINT64_MAX;
}

uint64_t SampleFromCdf(double cdf[], uint64_t size) {
  double r = UniformRandomDouble();
  uint64_t s = 0;
  for (uint64_t i = 0; i < size; i++) {
    s += (r > cdf[i]);
  }
  // TODO: use simd
  return s;
}

uint64_t SampleFromSurvival(double survival[], uint64_t size) {
  double r = UniformRandomDouble();
  uint64_t s = 0;
  for (uint64_t i = 0; i < size; i++) {
    s += (r < survival[i]);
  }
  return s;
}

template <typename OutputIterator>
uint64_t SampleFromPoisson(double lambda, OutputIterator begin,
                           OutputIterator end) {
  // requires lambda to be close to 1
  const uint64_t maxSample = 15;
  // calculate survival function
  double pdf[maxSample + 1];
  double survival[maxSample + 1];
  pdf[0] = exp(-lambda);
  for (uint64_t i = 1; i <= maxSample; i++) {
    pdf[i] = pdf[i - 1] * lambda / i;
  }
  survival[maxSample] =
      pdf[maxSample] * lambda / maxSample;  // estimation of tail
  for (uint64_t i = maxSample - 1; i > 0; i--) {
    survival[i] = survival[i + 1] + pdf[i + 1];
  }
  survival[0] = 1 - pdf[0];
  //   for (uint64_t i = 0; i <= maxSample; i++) {
  //     printf("%lu: %.10e\n", i, survival[i]);
  //   }
  uint64_t totalSample = 0;
  for (auto it = begin; it != end; it++) {
    totalSample += (*it = SampleFromSurvival(survival, maxSample + 1));
  }
  return totalSample;
}

uint64_t SampleFromBinomial(double p, uint64_t n) {
  const uint64_t maxSample = 15;
  double p_div_q = p / (1 - p);
  double pdf[maxSample + 1];
  double survival[maxSample + 1];

  pdf[0] = pow(1 - p, n);
  for (uint64_t i = 1; i <= maxSample; i++) {
    pdf[i] = pdf[i - 1] * (n - i + 1) * p_div_q / i;
  }

  survival[maxSample] = pdf[maxSample] * (n - maxSample) * p_div_q / maxSample;
  for (uint64_t i = maxSample - 1; i > 0; i--) {
    survival[i] = survival[i + 1] + pdf[i + 1];
  }
  survival[0] = 1 - pdf[0];
  //   for (uint64_t i = 0; i <= maxSample; i++) {
  //     printf("%lu: %.10e\n", i, survival[i]);
  //   }
  return SampleFromSurvival(survival, maxSample + 1);
}

template <typename OutputIterator>
uint64_t SampleFromBinomial(double p, uint64_t n, OutputIterator begin,
                            OutputIterator end) {
  uint64_t sum = 0;
  for (auto it = begin; it != end; it++) {
    sum += (*it = SampleFromBinomial(p, n));
  }
  return sum;
}

template <const bool approx = true>
struct NoReplaceSampler {
  uint64_t N, n, idx;
  static constexpr uint64_t cacheSize = approx ? 8 : 0;
  static_assert(cacheSize == 0 || (cacheSize & (cacheSize - 1)) == 0,
                "cacheSize must be power of 2");
  uint64_t cacheSamples[cacheSize];

  NoReplaceSampler(uint64_t N) : N(N), n(N), idx(0) {}

  NoReplaceSampler(uint64_t N, uint64_t n) : N(N), n(n), idx(0) {}

  uint64_t Sample() {
    if (n == 1) {
      --n;
      return N;
    }
    if constexpr (approx) {
      if (n > 1e5) {
        if (!(idx & (cacheSize - 1))) {
          uint64_t s = SampleFromPoisson((double)N / n, cacheSamples,
                                         cacheSamples + cacheSize);
          N -= s;
        }
        --n;
        ++idx;
        return cacheSamples[n & (cacheSize - 1)];
      } else {
        uint64_t s = SampleFromBinomial((double)1.0 / n, N);
        --n;
        N -= s;
        return s;
      }
    } else {
      uint64_t s = SampleFromBinomial((double)1.0 / n, N);
      --n;
      N -= s;
      return s;
    }
  }
};