#pragma once
#include <cfloat>
#include <cmath>
#include <vector>

#include "element.hpp"

/// This file contains some utility functions to help calculate parameters given
/// a desired success probability.

namespace Algorithm {

/// @brief add two numbers in log space
static double addLogs(double logA, double logB) {
  double bigger = std::max(logA, logB);
  double smaller = std::min(logA, logB);
  if (std::isinf(bigger)) {
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
    return logCkn + (double)k * log2(p) + (double)(n - k) * log2(1 - p);
  }
  return logCkn + (double)k * log2(p) -
         (double)(n - k) / log(2) * (p + p * p / 2 + p * p * p / 3);
}

/// @brief log of the binomial survival function
static double binomLogSf(uint64_t k, uint64_t n, double p) {
  double sf = -INFINITY;
  double pmf = binomLogPmf(k, n, p);
  double eps = pmf - 40;
  while (pmf > eps && k < n) {
    ++k;
    pmf += log2((double)(n - k + 1) / (double)k) + log2(p / (1 - p));
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
}  // namespace Algorithm