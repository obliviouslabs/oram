#pragma once
#include <vector>

#include "extemvector.hpp"

/// @brief Cache the first k entry of the array in internal memory, and the rest
/// in external memory encrypted. When accesses are made to the external
// memory data, cache the page in internal memory using a direct map cache.
namespace EM::CacheFrontVector {
using EM::ExtVector::EncryptType;
template <typename T,
          uint64_t page_size = std::max((1UL << 14) - 32, sizeof(T)),
          const EncryptType enc_type = EncryptType::ENCRYPT_AND_AUTH,
          const uint64_t ext_cache_bytes = (1UL << 18)>
struct Vector {
 private:
  // Set to log base 2, so compiler can optimize division in the direct map
  // cache
  static constexpr uint64_t cache_size =
      (1UL << GetLogBaseTwoConstexpr(ext_cache_bytes / page_size));
  using ExtVec = EM::ExtVector::Vector<T, page_size, enc_type, cache_size>;
  using IntVec = std::vector<T>;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  IntVec intVec;
  ExtVec* extVec = NULL;
  size_t cacheSize = 0;
  static uint64_t transformCacheBytesToCacheSize(uint64_t N,
                                                 uint64_t cacheBytes) {
    if (cacheBytes < ExtVec::GetMemoryUsage()) {
      throw std::runtime_error("Cache size too small");
    }
    return std::min(N, (cacheBytes - ExtVec::GetMemoryUsage()) / sizeof(T));
  }

 public:
  Vector() : intVec(0), cacheSize(0) {}

  explicit Vector(uint64_t N) : intVec(N), cacheSize(N) {}

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

  void SetSizeByCacheBytes(uint64_t N, uint64_t cacheBytes) {
    if (this->cacheSize != 0 || extVec != NULL) {
      throw std::runtime_error("SetSize can only be called on empty vector");
    }
    this->cacheSize = transformCacheBytesToCacheSize(N, cacheBytes);
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

  void SetSizeByCacheBytes(uint64_t N, uint64_t cacheBytes, const T& initVal) {
    if (this->cacheSize != 0 || extVec != NULL) {
      throw std::runtime_error("SetSize can only be called on empty vector");
    }
    this->cacheSize = transformCacheBytesToCacheSize(N, cacheBytes);
    intVec.resize(this->cacheSize, initVal);
    if (this->cacheSize < N) {
      extVec = new ExtVec(N - this->cacheSize, initVal);
    }
  }

  void InitDefault(const T& defaultVal) {
    for (size_t i = 0; i < cacheSize; ++i) {
      intVec[i] = defaultVal;
    }
    if (extVec)
      for (size_t i = 0; i < extVec->size(); ++i) {
        (*extVec)[i] = defaultVal;
      }
  }

  static uint64_t GetMinMemoryUsage() { return ExtVec::GetMemoryUsage(); }

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
      return extVec->Get(i - cacheSize);
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

  struct BatchAccessor {
    IntVec& intVec;
    size_t& cacheSize;
    typename ExtVec::BatchAccessor extAccessor;

    BatchAccessor(Vector& vec)
        : intVec(vec.intVec),
          cacheSize(vec.cacheSize),
          extAccessor(*(vec.extVec)) {
      Assert(vec.extVec);
    }

    uint32_t Prefetch(const uint64_t idx) {
      if (idx < cacheSize) {
        return -1;
      } else {
        return extAccessor.Prefetch(idx - cacheSize);
      }
    }

    void FlushRead() { extAccessor.FlushRead(); }

    T& At(uint32_t prefetchReceipt, uint64_t idx) {
      if (idx < cacheSize) {
        return intVec[idx];
      } else {
        return extAccessor.At(prefetchReceipt, idx - cacheSize);
      }
    }

    void FlushWrite() { extAccessor.FlushWrite(); }
  };
};
}  // namespace EM::CacheFrontVector