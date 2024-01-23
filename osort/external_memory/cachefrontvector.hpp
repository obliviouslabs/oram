#pragma once
#include <vector>

#include "extemvector.hpp"

// cache the first k entry of the arry in internal memory, and the rest in
// external memory (disk) encrypted. When accesses are made to the external
// memory cache the page in internal memory.
namespace EM::CacheFrontVector {
template <typename T,
          uint64_t page_size = std::max((1UL << 14) - 32, sizeof(T)),
          bool ENCRYPTED = true, bool AUTH = true>
struct Vector {
  static constexpr uint64_t item_per_page = page_size / sizeof(T);
  constexpr static bool useStdCopy = true;
  using ExtVec = EM::ExtVector::Vector<T, page_size, ENCRYPTED, AUTH, 1UL>;
  using IntVec = std::vector<T>;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  IntVec intVec;
  ExtVec extVec;
  size_t cacheSize;
  Vector(uint64_t N) : intVec(N), extVec(0), cacheSize(N) {}

  Vector(uint64_t N, uint64_t cacheSize)
      : intVec(std::min(cacheSize, N)),
        extVec(N - std::min(cacheSize, N)),
        cacheSize(std::min(cacheSize, N)) {}

  struct Iterator {
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using vector_type = Vector;
    constexpr static bool random_access = true;

    // Iterator constructors here...
    explicit Iterator(uint64_t ptr, Vector& vec) : m_ptr(ptr), vec_ptr(&vec) {}

    Iterator() : m_ptr(0), vec_ptr(NULL) {}

    reference operator*() {
      if (m_ptr < vec_ptr->cacheSize) {
        return vec_ptr->intVec[m_ptr];
      } else {
        return vec_ptr->extVec[m_ptr - vec_ptr->cacheSize];
      }
    }

    const_reference operator*() const {
      if (m_ptr < vec_ptr->cacheSize) {
        return vec_ptr->intVec[m_ptr];
      } else {
        return vec_ptr->extVec[m_ptr - vec_ptr->cacheSize];
      }
    }

    T* operator->() { return &(**this); }

    const T* operator->() const { return &(**this); }

    // Prefix increment
    Iterator& operator++() {
      ++m_ptr;
      return *this;
    }

    // Prefix decrement
    Iterator& operator--() {
      --m_ptr;
      return *this;
    }

    // Prefix increment
    Iterator& operator+=(int n) {
      m_ptr += n;
      return *this;
    }

    // Postfix increment
    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return a.m_ptr == b.m_ptr;
    };
    friend bool operator!=(const Iterator& a, const Iterator& b) {
      return a.m_ptr != b.m_ptr;
    };
    friend bool operator<(const Iterator& a, const Iterator& b) {
      return a.m_ptr < b.m_ptr;
    };
    friend bool operator<=(const Iterator& a, const Iterator& b) {
      return a.m_ptr <= b.m_ptr;
    };
    friend size_t operator-(const Iterator& it1, const Iterator& it2) {
      return it1.m_ptr - it2.m_ptr;
    }
    friend Iterator operator+(const Iterator& it, size_t size) {
      return Iterator(it.m_ptr + size, *it.vec_ptr);
    }
    friend Iterator operator-(const Iterator& it, size_t size) {
      return Iterator(it.m_ptr - size, *it.vec_ptr);
    }

   private:
    uint64_t m_ptr;
    Vector* vec_ptr;
  };

  Iterator begin() { return Iterator(0, *this); }

  Iterator end() { return Iterator(size(), *this); }

  size_t size() const { return intVec.size() + extVec.size(); }

  T& operator[](size_t i) {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      return extVec[i - cacheSize];
    }
  }

  const T& Get(size_t i) {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      return extVec.Get(i - cacheSize);
    }
  }
};
}  // namespace EM::CacheFrontVector