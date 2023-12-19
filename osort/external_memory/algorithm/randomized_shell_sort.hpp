#pragma once
#include <vector>

#include "common/utils.hpp"
#include "sort_def.hpp"
namespace EM::Algorithm {

class ShellSort {
 public:
  static const int C = 4;

  template <typename Iterator>
  static void compareExchange(Iterator i, Iterator j) {
    bool swapCond = (i < j) != (*i < *j);
    condSwap(swapCond, *i, *j);
  }

  template <typename Iterator>
  static void compareRegions(Iterator s, Iterator t, int offset) {
    std::vector<int> mate(offset);
    for (int count = 0; count < C; count++) {
      for (int i = 0; i < offset; i++) {
        mate[i] = i;
      }
      fisherYatesShuffle(mate.begin(), mate.end());
      for (int i = 0; i < offset; i++) {
        compareExchange(s + i, t + mate[i]);
      }
    }
  }

  template <typename Iterator>
  static void randomizedShellSort(Iterator begin, Iterator end) {
    int n = end - begin;
    for (int offset = n / 2; offset > 0; offset /= 2) {
      for (int i = 0; i < n - offset; i += offset) {
        compareRegions(begin + i, begin + i + offset, offset);
      }
      for (int i = n - offset; i >= offset; i -= offset) {
        compareRegions(begin + i - offset, begin + i, offset);
      }
      for (int i = 0; i < n - 3 * offset; i += offset) {
        compareRegions(begin + i, begin + i + 3 * offset, offset);
      }
      for (int i = 0; i < n - 2 * offset; i += offset) {
        compareRegions(begin + i, begin + i + 2 * offset, offset);
      }
      for (int i = 0; i < n; i += 2 * offset) {
        compareRegions(begin + i, begin + i + offset, offset);
      }
      for (int i = offset; i < n - offset; i += 2 * offset) {
        compareRegions(begin + i, begin + i + offset, offset);
      }
    }
  }
};

template <typename Iterator>
void RandomizedShellSort(Iterator begin, Iterator end) {
  ShellSort::randomizedShellSort(begin, end);
}

template <typename Vec>
void RandomizedShellSort(Vec& vec) {
  RandomizedShellSort(vec.begin(), vec.end());
}
}  // namespace EM::Algorithm