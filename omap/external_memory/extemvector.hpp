
#pragma once

#include <vector>

#include "external_memory/server/serverFrontend.hpp"
/// @brief External memory vector using a direct map / lru cache to swap pages.
namespace EM::ExtVector {
using EM::MemoryServer::EncryptType;
template <typename T,
          uint64_t page_size = std::max((1UL << 12) - 32, sizeof(T)),  // 16kB
          const EncryptType enc_type = EncryptType::ENCRYPT_AND_AUTH,
          uint64_t cache_size = SERVER__CACHE_SIZE>
struct Vector {
  static constexpr uint64_t item_per_page = page_size / sizeof(T);
  static constexpr bool AUTH = enc_type >= EncryptType::ENCRYPT_AND_AUTH;
  static constexpr bool ENCRYPTED = enc_type >= EncryptType::ENCRYPT;
  struct Page {
    T pages[item_per_page];
    using Encrypted_t = std::conditional_t<
        AUTH, FreshEncrypted<Page>,
        std::conditional_t<ENCRYPTED, Encrypted<Page>, NonEncrypted<Page>>>;
  };

  using Server = EM::MemoryServer::ServerFrontendInstance<
      Page, ::EM::Backend::MemServerBackend, enc_type, cache_size>;

  uint64_t N;
  Server server;

  static uint64_t GetMemoryUsage() { return Server::GetMemoryUsage(); }

  // https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
  struct Iterator {
    // Iterator tags here...
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = uint64_t;
    using page_idx_type = uint64_t;
    using page_offset_type = uint64_t;
    using reference = T&;
    using const_reference = const T&;
    using vector_type = Vector;

    // Iterator constructors here...
    explicit Iterator(pointer ptr, Vector& vec) : m_ptr(ptr), vec_ptr(&vec) {}

    Iterator() : m_ptr(0), vec_ptr(NULL) {}

    reference operator*() {
      Assert(m_ptr < vec_ptr->N);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / item_per_page;
      const size_t pageOffset = realIdx % item_per_page;
      return vec_ptr->server.Access(pageIdx).pages[pageOffset];
    }

    const_reference operator*() const {
      Assert(m_ptr < vec_ptr->N);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / item_per_page;
      const size_t pageOffset = realIdx % item_per_page;
      return vec_ptr->server.AccessReadOnly(pageIdx).pages[pageOffset];
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
    Vector* vec_ptr;
    pointer m_ptr;
    // Server* server;
  };

  // default:
  explicit Vector(uint64_t N_ = 0, typename Server::BackendType& _backend =
                                       *Backend::g_DefaultBackend)
      : N(N_), server(_backend, N_ / item_per_page + 1, makeDefaultPage()) {}
  Page makeDefaultPage() {
    T tdummy;
    return makeDefaultPage(tdummy);
  }
  Page makeDefaultPage(const T& defaultVal) {
    Page defaultPage;
    std::fill_n(defaultPage.pages, item_per_page, defaultVal);
    return defaultPage;
  }
  Vector(uint64_t N_, const T& defaultVal,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(N_),
        server(_backend, N_ / item_per_page + 1, makeDefaultPage(defaultVal)) {}
  // UNDONE: range and copy.

  Vector(Iterator begin, Iterator end,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(end - begin),
        server(_backend, (end - begin) / item_per_page + 1, makeDefaultPage()) {
    auto outIt = Iterator(0, *this);
    std::copy(begin, end, outIt);
  }

  Vector(Vector&& other) : N(other.N), server(std::move(other.server)) {}

  T& At(uint64_t index) {
    const size_t realIdx = index;
    const size_t pageIdx = realIdx / item_per_page;
    const size_t pageOffset = realIdx % item_per_page;
    return server.Access(pageIdx).pages[pageOffset];
  }

  const T& Get(uint64_t index) {
    const size_t realIdx = index;
    const size_t pageIdx = realIdx / item_per_page;
    const size_t pageOffset = realIdx % item_per_page;
    return server.AccessReadOnly(pageIdx).pages[pageOffset];
  }

  T& operator[](uint64_t index) { return At(index); }

  uint64_t size() const { return N; }

  Iterator begin() { return Iterator(0, *this); }

  Iterator end() { return Iterator(N, *this); }

  struct Reader {
    Iterator it;
    Iterator end;
    Reader(Iterator _begin, Iterator _end) : it(_begin), end(_end) {}

    const T& get() {
      const auto& it_const = it;
      return *it_const;
    }

    const T& read() {
      const T& val = get();
      ++it;
      return val;
    }
    bool eof() { return end <= it; }
  };

  struct Writer {
    Iterator it;
    Iterator end;
    Writer() {}
    Writer(Iterator _begin, Iterator _end) { init(_begin, _end); }
    void init(Iterator _begin, Iterator _end) {
      this->it = _begin;
      this->end = _end;
    }
    void write(const T& element) {
      *it = element;
      ++it;
    }
    bool eof() { return end <= it; }
    void flush() {}
  };

  struct BatchAccessor {
    typename Server::BatchAccessor serverAccessor;

    BatchAccessor(Vector& vec) : serverAccessor(vec.server, 64) {}

    uint32_t Prefetch(const uint64_t idx) {
      return serverAccessor.Prefetch(idx / item_per_page);
    }

    void FlushRead() { serverAccessor.FlushRead(); }

    T& At(uint32_t prefetchReceipt, uint64_t idx) {
      return serverAccessor.Access(prefetchReceipt).pages[idx % item_per_page];
    }

    void FlushWrite() { serverAccessor.FlushWrite(); }
  };
};

template <class InputIterator, class OutputIterator>
static OutputIterator Copy(InputIterator begin, InputIterator end,
                           OutputIterator to) {
  for (auto it = begin; it != end; ++it, ++to) {
    const auto& it_const = it;
    *to = *it_const;
  }
  return to;
}

}  // namespace EM::ExtVector