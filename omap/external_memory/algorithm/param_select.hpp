#pragma once
#include <cmath>
#include <vector>

#include "element.hpp"

/// This file contains some utility functions to help calculate parameters given
/// a desired success probability.

namespace EM::Algorithm {

/// @brief add two numbers in log space
static double addLogs(double logA, double logB) {
  double bigger = std::max(logA, logB);
  double smaller = std::min(logA, logB);
  if (bigger == -INFINITY) {
    return bigger;
  }
  if (bigger - smaller > 100) {
    return bigger;
  }
  return bigger + log2(1.0 + pow(2, (smaller - bigger)));
}

/// @brief log of the binomial coefficient
static double logCombin(uint64_t k, uint64_t n) {
  Assert(k <= n);
  if (k * 2 > n) {
    k = n - k;
  }

  double res = (lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1)) / log(2);
  return res;
}

/// @brief log of the binomial pmf
static double binomLogPmf(uint64_t k, uint64_t n, double p) {
  double logCkn = logCombin(k, n);

  if (p > 1e-5) {
    return logCkn + k * log2(p) + (n - k) * log2(1 - p);
  }
  return logCkn + k * log2(p) -
         (n - k) / log(2) * (p + p * p / 2 + p * p * p / 3);
}

/// @brief log of the binomial survival function
static double binomLogSf(uint64_t k, uint64_t n, double p) {
  double sf = -INFINITY;
  double pmf = binomLogPmf(k, n, p);
  double eps = pmf - 40;
  while (pmf > eps && k < n) {
    ++k;
    pmf += log2((double)(n - k + 1) / k) + log2(p / (1 - p));
    sf = addLogs(sf, pmf);
  }
  return sf;
}

/// @brief log of the binomial cdf
static double binomLogCdf(uint64_t k, uint64_t n, double p) {
  double sf = -INFINITY;
  double pmf = binomLogPmf(k, n, p);
  double eps = pmf - 40;
  while (pmf > eps) {
    sf = addLogs(sf, pmf);
    if (k == 0) {
      break;
    }
    pmf -= log2((double)(n - k + 1) / k) + log2(p / (1 - p));
    --k;
  }
  return sf;
}

/// @brief log of the hypergeometric pmf
static double hypergeomLogPmf(uint64_t k, uint64_t M, uint64_t n, uint64_t N) {
  Assert(N <= M);
  if (k > N || k > n) {
    return -INFINITY;
  }

  double logCkn = logCombin(k, n);
  double logCN_kM_n = logCombin(N - k, M - n);
  static size_t cachedM = 0;
  static size_t cachedN = 0;
  static size_t cachedLogMN = 0;
  double logCMN;
  if (M == cachedM && N == cachedN) {
    logCMN = cachedLogMN;
  } else {
    logCMN = logCombin(N, M);
    cachedM = M;
    cachedN = N;
    cachedLogMN = logCMN;
  }
  return logCkn + logCN_kM_n - logCMN;
}

/// @brief log of the hypergeometric survival function
static double hypergeomLogSf(uint64_t k, uint64_t M, uint64_t n, uint64_t N) {
  double sf = -INFINITY;
  double pmf = hypergeomLogPmf(k, M, n, N);
  double eps = pmf - 40;
  while (pmf > eps && k < n) {
    ++k;
    pmf += log2((double)(n - k + 1) / k) -
           log2((double)(M - n - N + k) / (N - k + 1));
    sf = addLogs(sf, pmf);
  }
  return sf;
}

/// @brief calculate lower bound
/// @tparam T type of the value
/// @tparam type of the lambda function
/// @param left left bound
/// @param right right bound
/// @param satisfy lambda function that checks the condition
/// @param prec precision
/// @return the smallest value in the range [left, right) that satisfies the
/// condition
template <typename T, class Check>
static T lowerBound(T left, T right, const Check& satisfy, T prec = 1) {
  while (true) {
    T mid = (left + right) / 2;
    if (satisfy(mid)) {
      right = mid;
      if (right - left <= prec) {
        return mid;
      }
    } else {
      left = mid;
      if (right - left <= prec) {
        return right;
      }
    }
  }
}
}  // namespace EM::Algorithm