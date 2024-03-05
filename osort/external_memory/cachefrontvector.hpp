#pragma once
#include <vector>

#include "extemvector.hpp"

// cache the first k entry of the arry in internal memory, and the rest in
// external memory (disk) encrypted. When accesses are made to the external
// memory cache the page in internal memory.
namespace EM::CacheFrontVector {
template <typename T,
          uint64_t page_size = std::max((1UL << 14) - 32, sizeof(T)),
          bool ENCRYPTED = true, bool AUTH = true,
          const uint64_t ext_cache_bytes = (1UL << 16)>
struct Vector {
  static constexpr uint64_t item_per_page = page_size / sizeof(T);
  constexpr static bool useStdCopy = true;
  using ExtVec =
      EM::ExtVector::Vector<T, page_size, ENCRYPTED, AUTH,
                            std::max(1UL, ext_cache_bytes / page_size)>;
  using IntVec = std::vector<T>;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  IntVec intVec;
  ExtVec* extVec = NULL;
  size_t cacheSize = 0;

  Vector() : intVec(0), cacheSize(0) {}

  Vector(uint64_t N) : intVec(N), cacheSize(N) {}

  Vector(uint64_t N, uint64_t cacheSize) { SetSize(N, cacheSize); }

  Vector(uint64_t N, uint64_t cacheSize, const T& initVal) {
    SetSize(N, cacheSize, initVal);
  }

  ~Vector() {
    if (extVec) {
      delete extVec;
      extVec = NULL;
    }
  }

  void SetSize(uint64_t N) {
    if (cacheSize != 0 || extVec != NULL) {
      throw std::runtime_error("SetSize can only be called on empty vector");
    }
    intVec.resize(N);
  }

  void SetSize(uint64_t N, uint64_t cacheSize) {
    if (this->cacheSize != 0 || extVec != NULL) {
      throw std::runtime_error("SetSize can only be called on empty vector");
    }
    this->cacheSize = std::min(cacheSize, N);
    intVec.resize(this->cacheSize);
    if (this->cacheSize < N) {
      extVec = new ExtVec(N - this->cacheSize);
    }
  }

  void SetSize(uint64_t N, uint64_t cacheSize, const T& initVal) {
    if (this->cacheSize != 0 || extVec != NULL) {
      throw std::runtime_error("SetSize can only be called on empty vector");
    }
    this->cacheSize = std::min(cacheSize, N);
    intVec.resize(this->cacheSize, initVal);
    if (this->cacheSize < N) {
      extVec = new ExtVec(N - this->cacheSize, initVal);
    }
  }

  static uint64_t GetMemoryUsage(uint64_t N, uint64_t cacheSize) {
    if (cacheSize < N) {
      return cacheSize * sizeof(T) + ExtVec::GetMemoryUsage();
    } else {
      return cacheSize * sizeof(T);
    }
  }

  uint64_t GetMemoryUsage() const {
    return Vector::GetMemoryUsage(size(), cacheSize);
  }

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
        Assert(vec_ptr->extVec);
        return (*vec_ptr->extVec)[m_ptr - vec_ptr->cacheSize];
      }
    }

    const_reference operator*() const {
      if (m_ptr < vec_ptr->cacheSize) {
        return vec_ptr->intVec[m_ptr];
      } else {
        Assert(vec_ptr->extVec);
        return (*vec_ptr->extVec)[m_ptr - vec_ptr->cacheSize];
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
    Iterator& operator+=(uint64_t n) {
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

  size_t size() const {
    if (extVec == NULL) {
      return cacheSize;
    } else {
      return cacheSize + extVec->size();
    }
  }

  T& operator[](size_t i) {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      Assert(extVec);
      return (*extVec)[i - cacheSize];
    }
  }

  const T& operator[](size_t i) const {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      Assert(extVec);
      return (*extVec)[i - cacheSize];
    }
  }

  const T& Get(size_t i) {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      Assert(extVec);
      return extVec->Get(i - cacheSize);
    }
  }

  T& GetMutable(size_t i) {
    if (i < cacheSize) {
      return intVec[i];
    } else {
      Assert(extVec);
      return extVec->At(i - cacheSize);
    }
  }

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

    INLINE const T& get() const {
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
}  // namespace EM::CacheFrontVector