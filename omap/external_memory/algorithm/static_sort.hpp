#pragma once
#include "sort_def.hpp"

/// This file contains a static sort class that uses a compile time generated.

// https://stackoverflow.com/questions/19790522/very-fast-sorting-of-fixed-length-arrays-using-comparator-networks
/// @brief A Functor class to create a sort for fixed sized arrays/containers
/// with a compile time generated Bose-Nelson sorting network.
/// @tparam NumElements  The number of elements in the array or container to
/// sort
/// @tparam T            The element type.
/// @tparam Compare      A comparator functor class that returns true if lhs <
/// rhs.
namespace EM::Algorithm {
template <unsigned NumElements, class Compare = void>
class StaticSort {
  template <class A, class C>
  struct Swap {
    inline Swap(A &a, const int &i0, const int &i1) {
      bool cond = Compare()(*(a + i0), *(a + i1));
      condSwap(cond, *(a + i0), *(a + i1));
    }
  };

  template <class A>
  struct Swap<A, void> {
    inline Swap(A &a, const int &i0, const int &i1) {
      bool cond = *(a.second + i0) > *(a.second + i1);
      condSwap(cond, *(a.first + i0), *(a.first + i1));
      condSwap<false>(cond, *(a.second + i0), *(a.second + i1));
    }
  };

  template <class A, class C, int I, int J, int X, int Y>
  struct PB {
    inline PB(A &a) {
      enum {
        L = X >> 1,
        M = (X & 1 ? Y : Y + 1) >> 1,
        IAddL = I + L,
        XSubL = X - L
      };
      PB<A, C, I, J, L, M> p0(a);
      PB<A, C, IAddL, J + M, XSubL, Y - M> p1(a);
      PB<A, C, IAddL, J, XSubL, M> p2(a);
    }
  };

  template <class A, class C, int I, int J>
  struct PB<A, C, I, J, 1, 1> {
    inline PB(A &a) { Swap<A, C> s(a, I - 1, J - 1); }
  };

  template <class A, class C, int I, int J>
  struct PB<A, C, I, J, 1, 2> {
    inline PB(A &a) {
      Swap<A, C> s0(a, I - 1, J);
      Swap<A, C> s1(a, I - 1, J - 1);
    }
  };

  template <class A, class C, int I, int J>
  struct PB<A, C, I, J, 2, 1> {
    inline PB(A &a) {
      Swap<A, C> s0(a, I - 1, J - 1);
      Swap<A, C> s1(a, I, J - 1);
    }
  };

  template <class A, class C, int I, int M, bool Stop = false>
  struct PS {
    inline PS(A &a) {
      enum { L = M >> 1, IAddL = I + L, MSubL = M - L };
      PS<A, C, I, L, (L <= 1)> ps0(a);
      PS<A, C, IAddL, MSubL, (MSubL <= 1)> ps1(a);
      PB<A, C, I, IAddL, L, MSubL> pb(a);
    }
  };

  template <class A, class C, int I, int M>
  struct PS<A, C, I, M, true> {
    inline PS(A &a) {}
  };

 public:
  /**
   * Sorts the array/container arr.
   * \param  arr  The array/container to be sorted.
   */
  template <class Container>
  inline void operator()(Container &arr) const {
    PS<Container, Compare, 1, NumElements, (NumElements <= 1)> ps(arr);
  };
};
}  // namespace EM::Algorithm