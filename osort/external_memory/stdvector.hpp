#pragma once
#include <vector>

#include "common/defs.hpp"
/**
 * A wrapper for std::vector with the same interface as the other containers
 * */
template <typename T>
struct StdVector : public std::vector<T> {
  constexpr static bool useStdCopy = true;
  constexpr static size_t item_per_page = 1;

  struct Iterator : public std::vector<T>::iterator {
    using vector_type = StdVector<T>;
    constexpr static bool random_access = true;

    Iterator() : std::vector<T>::iterator() {}

    // Constructor to initialize the base iterator and custom attribute
    Iterator(const typename std::vector<T>::iterator& it)
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
    size_t get_page_offset() { return 0; }
  };

  StdVector(uint64_t N) : std::vector<T>(N) {}
  StdVector(Iterator begin, Iterator end) : std::vector<T>(begin, end) {}

  StdVector(uint64_t N, const T& val) : std::vector<T>(N, val) {}

  Iterator begin() { return Iterator(std::vector<T>::begin()); }

  Iterator end() { return Iterator(std::vector<T>::end()); }

  INLINE size_t size() { return std::vector<T>::size(); }

  struct Reader {
    using value_type = T;
    using iterator_type = Iterator;
    Iterator it;
    Iterator end;

    Reader() {}

    Reader(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
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

  struct PrefetchReader : public Reader {
    using value_type = T;
    using iterator_type = Iterator;
    PrefetchReader() {}

    PrefetchReader(Iterator _begin, Iterator _end, uint32_t _auth = 0,
                   uint64_t _heapSize = DEFAULT_HEAP_SIZE)
        : Reader(_begin, _end, _auth) {}
  };

  struct Writer {
    using value_type = T;
    using iterator_type = Iterator;
    Iterator it;
    Iterator end;
    Writer() {}

    Writer(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
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