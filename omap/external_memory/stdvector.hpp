#pragma once
#include <vector>

#include "common/defs.hpp"
/**
 * A wrapper for std::vector with the same interface as the other vector
 * containers
 * */
template <typename T>
struct StdVector : public std::vector<T> {
  struct Iterator : public std::vector<T>::iterator {
    using vector_type = StdVector<T>;
    using iterator_category = std::random_access_iterator_tag;

    Iterator() : std::vector<T>::iterator() {}

    // Constructor to initialize the base iterator and custom attribute
    explicit Iterator(const typename std::vector<T>::iterator& it)
        : std::vector<T>::iterator(it) {}

    friend Iterator operator+(const Iterator& it, size_t size) {
      return Iterator(typename std::vector<T>::iterator(it) + size);
    }

    friend Iterator operator-(const Iterator& it, size_t size) {
      return Iterator(typename std::vector<T>::iterator(it) - size);
    }

    friend Iterator operator+(const Iterator& it, int64_t size) {
      return Iterator(typename std::vector<T>::iterator(it) + size);
    }

    friend Iterator operator-(const Iterator& it, int64_t size) {
      return Iterator(typename std::vector<T>::iterator(it) - size);
    }

    friend Iterator operator+(const Iterator& it, int32_t size) {
      return Iterator(typename std::vector<T>::iterator(it) + size);
    }

    friend Iterator operator-(const Iterator& it, int32_t size) {
      return Iterator(typename std::vector<T>::iterator(it) - size);
    }

    friend Iterator operator+(const Iterator& it, uint32_t size) {
      return Iterator(typename std::vector<T>::iterator(it) + size);
    }

    friend Iterator operator-(const Iterator& it, uint32_t size) {
      return Iterator(typename std::vector<T>::iterator(it) - size);
    }
  };
  StdVector() : std::vector<T>() {}
  explicit StdVector(uint64_t N) : std::vector<T>(N) {}
  StdVector(Iterator begin, Iterator end) : std::vector<T>(begin, end) {}

  StdVector(uint64_t N, const T& val) : std::vector<T>(N, val) {}
  StdVector(const StdVector& other) : std::vector<T>(other) {}
  StdVector(StdVector&& other) : std::vector<T>(std::move(other)) {}

  Iterator begin() { return Iterator(std::vector<T>::begin()); }

  Iterator end() { return Iterator(std::vector<T>::end()); }

  void resize(uint64_t N) { std::vector<T>::resize(N); }
  void SetSize(uint64_t N, uint64_t _ignore1 = 0, uint64_t _ignore2 = 0) {
    resize(N);
  }
  INLINE size_t size() { return std::vector<T>::size(); }

  struct Reader {
    using value_type = T;
    using iterator_type = Iterator;
    Iterator it;
    Iterator end;

    Reader() {}

    Reader(Iterator _begin, Iterator _end) : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end) {
      it = _begin;
      end = _end;
    }

    INLINE T& get() {
      Assert(!eof());
      return *it;
    }

    INLINE const T& read() {
      const T& val = get();
      ++it;
      return val;
    }

    INLINE bool eof() { return end <= it; }

    size_t size() { return end - it; }
  };

  struct Writer {
    using value_type = T;
    using iterator_type = Iterator;
    Iterator it;
    Iterator end;
    Writer() {}

    Writer(Iterator _begin, Iterator _end) : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end) {
      it = _begin;
      end = _end;
    }

    INLINE void write(const T& element) {
      *it = element;
      ++it;
    }

    INLINE bool eof() { return end <= it; }

    size_t size() { return end - it; }

    INLINE void flush() {}
  };
};