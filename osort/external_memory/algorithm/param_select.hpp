#pragma once
#include <cmath>
#include <vector>

#include "sort_def.hpp"

/// This file contains the parameter solver for flex-way algorithms.

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

/// @brief calculate failure probability of flex-way butterfly sort
/// @param Z bucket size
/// @param N number of real elements
/// @param numTotalBucket number of total buckets
/// @return log of the failure probability
static double failureProbButterflySort(size_t Z, size_t N,
                                       size_t numTotalBucket) {
  size_t NTotal = numTotalBucket * Z;
  double failProb = -INFINITY;
  size_t numRealPerBucket = (N - 1) / numTotalBucket + 1;
  for (size_t numBucket = numTotalBucket; numBucket >= 2; numBucket /= 2) {
    size_t n = numRealPerBucket * numBucket;  // num real on first level
    if (Z >= n) {
      break;
    }
    double p = 1.0 / numBucket;
    double bucketFailProb = binomLogSf(Z, n, p);
    double layerFailProb = log2(numTotalBucket) + bucketFailProb;
    failProb = addLogs(failProb, layerFailProb);
  }
  return failProb;
}

/// @brief judge whether the failure probability of flex-way butterfly sort is
/// within the target
/// @param Z bucket size
/// @param N number of real elements
/// @param B granularity of preshuffling
/// @param numTotalBucket number of total buckets
/// @param pMax maximum number of partition
/// @param sampleRatio sampling rate
/// @param target target failure probability in log scale
/// @return whether the failure probability is within the target
static bool failureProbKWayDistriSortWithinTarget(size_t Z, size_t N, size_t B,
                                                  size_t numTotalBucket,
                                                  size_t pMax,
                                                  double sampleRatio,
                                                  double target = -60) {
  size_t NTotal = numTotalBucket * Z;
  double failProb = -INFINITY;
  size_t numRealPerBucket = (N - 1) / numTotalBucket + 1;
  size_t layer = log2(numTotalBucket);
  for (size_t numBucket = pMax; numBucket >= 2; numBucket /= 2) {
    size_t n = numRealPerBucket * numBucket;  // num real on first level
    if (Z >= n) {
      break;
    }
    auto not_satisfy_start = [&](uint64_t t1) {
      double probOverflow = hypergeomLogSf(Z / B, N / B, t1 / B, n / B);
      double probAnyOverflow = probOverflow + log2(numBucket);
      return probAnyOverflow >= target - layer - 3;
    };
    size_t t1_start =
        lowerBound(N / numBucket, NTotal / numBucket, not_satisfy_start) - 1;

    auto satisfy_end = [&](int64_t t1) {
      double probMoreThanT1BtwPivots =
          binomLogCdf(ceil(N * sampleRatio / numBucket), t1 - 1, sampleRatio) +
          log2(N);
      // probability that the number of elements between any two
      // pivots are more than t1
      double probMoreThanT1BtwAnyPivots =
          probMoreThanT1BtwPivots + log2(numBucket);
      return probMoreThanT1BtwAnyPivots < target - layer - 3;
    };
    size_t t1_end = lowerBound(N / numBucket, NTotal / numBucket, satisfy_end);
    if (t1_start >= t1_end) {
      t1_start = t1_end = (t1_start + t1_end) / 2;
    }

    // add failure prob on two ends
    // either end < target -layer - 3
    double qLayer = target - layer - 2;

    size_t step = std::max(1UL, (t1_end - t1_start) / 100);
    for (size_t t1 = t1_start; t1 < t1_end; t1 += step) {
      // P[max(element between any two pivots) = t1]
      // <= P[there exist two pivots sandwiching t1 elements]
      // <= P[there exists a consecutive chunk of t1-1 elements where
      // total*sampleratio/numbucket elements are sampled and both ends are
      // sampled]
      // <= 2^prob
      double qi = binomLogPmf(ceil(N * sampleRatio / numBucket) - 1, t1 - 1,
                              sampleRatio) +
                  log2(sampleRatio);
      double qj = hypergeomLogSf(Z / B, N / B, (t1 + step) / B, n / B);
      qLayer = addLogs(qLayer, qi + qj + log2(numBucket) +
                                   log2(std::min(step, t1_end - t1)));
    }
    // qLayer += log2(log2(std::min(way, p)));
    failProb = addLogs(failProb, qLayer);
    if (failProb > target) {
      return false;
    }
    double ub = addLogs(failProb, qLayer + log2(numBucket / 2));
    if (ub < target) {
      return true;
    }
  }

  return failProb < target;
}

/// @brief solve the optimal way of flex-way butterfly/distribution sort
class ButterflyWaySolver {
 private:
  // options of mergesplit ways
  std::vector<uint64_t> options;
  // amortized cost on each bucket in option-way mergesplit
  std::vector<double> unitCosts;
  // amortized cost to read&write a bucket from/to ext mem
  double ioCost;
  // maximum number of ways that can be done in internal memory
  uint64_t maxWayInternal;

  double solverHelper(std::vector<uint64_t>& optimalWays, size_t optionCount,
                      uint64_t targetProduct, uint64_t remainWayInternal) {
    optimalWays.clear();
    if (targetProduct <= 1) {
      return 0;
    }
    if (!optionCount) {
      return INFINITY;
    }
    std::vector<uint64_t> childrenOptimalWays;
    double minCost = INFINITY;
    for (size_t optionIdx = 0; optionIdx < optionCount; ++optionIdx) {
      uint64_t option = options[optionIdx];
      if (option <= remainWayInternal) {
        // for efficiency only allow to use options <= current option
        double costChildren = solverHelper(childrenOptimalWays, optionIdx + 1,
                                           divRoundUp(targetProduct, option),
                                           remainWayInternal / option);
        uint64_t product = option;
        for (uint64_t way : childrenOptimalWays) {
          product *= way;
        }
        double cost = costChildren * option + unitCosts[optionIdx] * product;
        if (cost < minCost) {
          minCost = cost;
          optimalWays.resize(childrenOptimalWays.size() + 1);
          std::copy(childrenOptimalWays.begin(), childrenOptimalWays.end(),
                    optimalWays.begin());
          optimalWays.back() = option;
        }
        if (option == targetProduct || option == maxWayInternal) {
          return minCost;  // no need to try next option
        }
      }
    }
    if (options[optionCount - 1] > remainWayInternal) {
      // start new layer, reset all options and remainWayInternal
      double costWithoutIO = solverHelper(childrenOptimalWays, options.size(),
                                          targetProduct, maxWayInternal);
      uint64_t product = 1;
      for (uint64_t way : childrenOptimalWays) {
        product *= way;
      }
      double cost = costWithoutIO + product * ioCost;
      if (cost < minCost) {
        minCost = cost;
        optimalWays.resize(childrenOptimalWays.size());
        std::copy(childrenOptimalWays.begin(), childrenOptimalWays.end(),
                  optimalWays.begin());
      }
    }
    return minCost;
  }

 public:
  ButterflyWaySolver(const std::vector<uint64_t>& options,
                     const std::vector<double>& unitCosts, double ioCost,
                     uint64_t maxWayInternal)
      : options(options),
        unitCosts(unitCosts),
        ioCost(ioCost),
        maxWayInternal(maxWayInternal) {
    Assert(unitCosts.size() == options.size());
    std::sort(this->options.begin(), this->options.end());
  }

  double solve(std::vector<std::vector<uint64_t>>& output,
               uint64_t targetProduct) {
    std::vector<uint64_t> optimalWays;
    double cost = solverHelper(optimalWays, options.size(), targetProduct,
                               maxWayInternal);
    output.clear();
    uint64_t accWays = 1;
    auto begin = optimalWays.begin();
    auto it = begin;
    for (; it != optimalWays.end(); ++it) {
      if (accWays * (*it) > maxWayInternal) {
        if (begin != it) {
          output.emplace_back(begin, it);
        }
        begin = it;
        accWays = *it;
      } else {
        accWays *= (*it);
      }
    }
    output.emplace_back(begin, it);
    return cost;
  }
};

struct KWayButterflyParams {
  std::vector<std::vector<uint64_t>> ways;  // ways of MergeSplit on each layer
  size_t Z;                                 // bucket size
  size_t totalBucket;                       // total number of buckets
};

struct DistriParams {
  double samplingRatio;                     // sampling rate
  double slackSampling;                     // slack of sampling in each batch
  std::vector<std::vector<uint64_t>> ways;  // ways of MergeSplit on each layer
  size_t Z;                                 // bucket size
  size_t totalPartition;                    // product of all ways
  size_t totalBucket;  // total number of buckets, should be >= totalPartition
  size_t minTotalBucket;  // minimum total number of buckets to satisfy failure
                          // probability
};

/// @brief find the optimal parameters for flex-way butterfly o-sort
/// @param N number of elements
/// @param M number of elements that fit in internal memory
/// @param target failure probability target in log scale
/// @return optimal parameters
static KWayButterflyParams bestKWayButterflyParams(size_t N,
                                                   size_t M = 0x7000000 / 136,
                                                   int64_t b = 128,
                                                   double target = -60) {
  double minCost = INFINITY;
  KWayButterflyParams optimalParams;
  const std::vector<uint64_t> choices = {2, 3, 4, 5, 6, 7, 8};
  const double bitonicCostFactor = 7.2e-6;

  const std::vector<std::vector<double>> costs = {
      {6.889e-02, 1.698e-07, 1.905e-07}, {5.630e-02, 9.017e-06, 1.781e-07},
      {6.714e-02, 1.122e-05, 1.691e-07}, {7.285e-02, 1.247e-05, 1.734e-07},
      {5.692e-02, 1.221e-05, 1.825e-07}, {5.723e-02, 1.303e-05, 1.878e-07},
      {5.684e-02, 1.345e-05, 1.892e-07}};

  std::vector<double> costZ(7);

#ifdef DISK_IO
  static const double ioCostPerElement = 0.00096857;
#else
  static const double ioCostPerElement = 0.00061467;
#endif
  for (int i = 0; i < 7; ++i) {
    size_t logZ = 8 + i;
    size_t Z = 1UL << logZ;
    size_t M0 = M - divRoundUp((b + 2) * Z * 8, b);
    if (Z > N / 2 || 16 * Z > M0) {
      if (i == 0) {
        printf("Input size / memory size too small\n");
        abort();
      }
      break;
    }
    auto satisfy = [=](size_t bucketCount) {
      return failureProbButterflySort(Z, N, bucketCount) < target;
    };
    size_t minBucketCount = lowerBound(N / Z + 1, 3 * N / Z, satisfy);
    size_t maxWayInternal = M0 / Z;
    for (int j = 0; j < 7; ++j) {
      costZ[j] = costs[j][0] + (costs[j][1] + b * costs[j][2]) * Z * logZ;
    }
    ButterflyWaySolver solver(choices, costZ, ioCostPerElement * Z,
                              maxWayInternal);
    KWayButterflyParams params;
    double cost = solver.solve(params.ways, minBucketCount);
    uint64_t actualBucketCount = 1;
    for (auto& vec : params.ways) {
      for (auto way : vec) {
        actualBucketCount *= way;
      }
    }
    cost += bitonicCostFactor * Z * logZ * logZ * actualBucketCount;
    if (cost < minCost) {
      minCost = cost;
      optimalParams = params;
      optimalParams.Z = Z;
      optimalParams.totalBucket = actualBucketCount;
    }
  }
  return optimalParams;
}

/// @brief find the optimal parameters for flex-way distri o-sort
/// @param N number of elements
/// @param M number of elements that fit in internal memory
/// @param B number of elements per page, set to 1 if already pre-shuffled
/// @param target failure probability target in log scale
/// @return optimal parameters
static DistriParams bestDistriParams(size_t N, size_t M = 0x7000000 / 136,
                                     size_t B = 1, int64_t b = 128,
                                     double target = -60) {
  double minCost = INFINITY;
  DistriParams optimalParams;
  const std::vector<uint64_t> choices = {2, 3, 4, 5, 6, 7, 8};
  const double bitonicCostFactor = 7.2e-6;
  // const std::vector<std::vector<double>> costs16 = {
  //     {0.027208, 0.062339, 0.074877, 0.083373, 0.084987, 0.090973, 0.094594},
  //     {0.053427, 0.131716, 0.155852, 0.164248, 0.113467, 0.120657, 0.126136},
  //     {0.067885, 0.178683, 0.209582, 0.230375, 0.233684, 0.250123, 0.262026},
  //     {0.148892, 0.374595, 0.436384, 0.477227, 0.482734, 0.514162, 0.535860},
  //     {0.304131, 0.773884, 0.901490, 0.984308, 0.999014, 1.059123, 1.103032},
  //     {0.624599, 1.621319, 1.873467, 2.033733, 2.064392, 2.182470, 2.274683},
  //     {1.301665, 3.375710, 3.884231, 4.207232, 4.271370, 4.504563, 4.711345},
  //     {3.247753, 7.079485, 8.103368, 8.672383, 8.987384, 9.347671, 9.674656}};
  // const std::vector<std::vector<double>> costs128 = {
  //     {0.040987, 0.072334, 0.077112, 0.083998, 0.085417, 0.090983, 0.093360},
  //     {0.088015, 0.155585, 0.163523, 0.177144, 0.180525, 0.192080, 0.197212},
  //     {0.200702, 0.329450, 0.343935, 0.379935, 0.394584, 0.422252, 0.437080},
  //     {0.413447, 0.730196, 0.780530, 0.855178, 0.875196, 0.917903, 0.946604},
  //     {0.885438, 1.634678, 1.715012, 1.831199, 1.864355, 1.961980, 1.999063},
  //     {2.571125, 3.490637, 3.658499, 3.926035, 4.011425, 4.294402, 4.435121},
  //     {5.754265, 7.454927, 7.822982, 8.327188, 8.439162, 8.859011, 8.953003},
  //     {12.228781, 15.865265, 16.463026, 17.766671, 17.935565, 18.941199,
  //      19.194945}};

  // std::vector<std::vector<double>> costs(
  //     costs128.size(), std::vector<double>(costs128[0].size()));
  // for (size_t row = 0; row < costs128.size(); ++row) {
  //   for (size_t col = 0; col < costs128[row].size(); ++col) {
  //     costs[row][col] =
  //         (costs128[row][col] * (b - 16) + costs16[row][col] * (128 - b)) /
  //         112 * 136 / (b + 8);
  //   }
  // }

  const std::vector<std::vector<double>> costs = {
      {1.000e-02, 5.000e-06, 1.791e-07}, {7.015e-02, 9.300e-06, 1.762e-07},
      {8.854e-02, 1.160e-05, 1.680e-07}, {1.137e-01, 1.236e-05, 1.801e-07},
      {7.767e-02, 1.308e-05, 1.807e-07}, {1.000e-01, 1.347e-05, 1.923e-07},
      {1.092e-01, 1.421e-05, 1.908e-07}};
  std::vector<double> costZ(7);
#ifdef DISK_IO
  static const double ioCostPerElement = 0.00096857;
#else
  static const double ioCostPerElement = 0.00061467;
#endif
  for (double samplingRatio = 0.01; samplingRatio < 0.031;
       samplingRatio += 0.005) {
    for (int i = 0; i < 12; ++i) {
      size_t logZ = 8 + i;
      size_t Z = 1UL << logZ;
      size_t M0 =
          M - divRoundUp((b + 2) * Z * 8, b) -
          4 * N / M;  // minus space used for mergesplit and storing pivots
      if (Z > N / 2 || 16 * Z > M0) {
        if (i == 0) {
          printf("Input size / memory size too small\n");
          abort();
        }
        break;
      }
      size_t maxFinalPartSize = M0 * b / (b + 4);
      auto satisfy = [=](size_t bucketCount) {
        return failureProbKWayDistriSortWithinTarget(
            Z, N, B, bucketCount, divRoundUp(N, maxFinalPartSize),
            samplingRatio, target);
      };

      size_t minBucketCount = lowerBound(N / Z + 1, 5 * N / Z, satisfy);
      size_t minTotalWay = divRoundUp(minBucketCount, maxFinalPartSize / Z);
      size_t maxWayInternal = M0 / Z;
      for (int j = 0; j < 7; ++j) {
        costZ[j] = costs[j][0] + (costs[j][1] + b * costs[j][2]) * Z * logZ;
      }
      ButterflyWaySolver solver(choices, costZ, ioCostPerElement * Z,
                                maxWayInternal);
      DistriParams params;
      double cost = solver.solve(params.ways, minTotalWay);
      uint64_t totalPartition = 1;
      for (auto& vec : params.ways) {
        for (auto way : vec) {
          totalPartition *= way;
        }
      }
      auto satisfy2 = [=](size_t bucketCountFinalLayer) {
        return failureProbKWayDistriSortWithinTarget(
            Z, N, B, bucketCountFinalLayer * totalPartition, totalPartition,
            samplingRatio, target);
      };
      size_t bucketCountFinalLayer =
          lowerBound(divRoundUp(minBucketCount, totalPartition),
                     maxFinalPartSize / Z, satisfy2);
      cost *= bucketCountFinalLayer;
      size_t partitionSize = bucketCountFinalLayer * Z;
      double logPartitionSize = log2(partitionSize);
      cost += bitonicCostFactor * partitionSize * logPartitionSize *
              logPartitionSize * totalPartition;
      auto samplingSatisify = [=](double slack) {
        return binomLogSf(slack * M * samplingRatio, M, samplingRatio) +
                   log2(divRoundUp(N, M)) <
               target - 4;
      };
      double slackSampling = lowerBound(1.0, 1.5, samplingSatisify, 1e-4);
      cost *= 1 + samplingRatio * slackSampling;
      if (cost < minCost) {
        minCost = cost;
        optimalParams = params;  // copy ways
        optimalParams.Z = Z;
        optimalParams.totalPartition = totalPartition;
        optimalParams.totalBucket = totalPartition * bucketCountFinalLayer;
        optimalParams.minTotalBucket = minBucketCount;
        optimalParams.samplingRatio = samplingRatio;
        optimalParams.slackSampling = slackSampling;
      }
    }
  }

  return optimalParams;
}
}  // namespace EM::Algorithm