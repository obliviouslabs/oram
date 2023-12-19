
#pragma once

#include <vector>

#include "common/encrypted.hpp"
#include "external_memory/server/serverFrontend_dynamic.hpp"

///@brief External memory vector that allows page size to be dynamically
/// decided. Only supports CopyIn/CopyOut.
namespace EM::DynamicPageVector {
template <typename T, bool ENCRYPTED = true, bool AUTH = true>
struct Vector {
  static constexpr bool encrypted = ENCRYPTED;  // whether to use encryption
  static constexpr bool auth = AUTH;  // whether to use authenticated encryption
  constexpr static bool useStdCopy = false;

  uint64_t referenceCount;
  uint64_t item_per_page;

  // using Server =
  // EM::MemoryServer::NonCachedServerFrontendInstance<Page,::EM::Backend::MemServerBackend,ENCRYPTED>;

  using Server = EM::MemoryServer::DynamicServerFrontendInstance<
      ::EM::Backend::MemServerBackend, ENCRYPTED, AUTH>;

  // fileserver here.
  //
  uint64_t N;
  Server server;

  // https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
  struct Iterator {
    // Iterator tags here...
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = uint64_t;
    using page_idx_type = uint64_t;
    using page_offset_type = uint64_t;
    using reference = T&;
    using const_reference = const T&;
    using vector_type = Vector;
    constexpr static bool random_access = false;
    // Iterator constructors here...
    explicit Iterator(pointer ptr, Vector& vec) : m_ptr(ptr), vec_ptr(&vec) {}

    Iterator() : m_ptr(0), vec_ptr(NULL) {}

    Iterator(pointer ptr) : m_ptr(ptr), vec_ptr(NULL) {}

    // should not be called unless absolutely necessary
    T operator*() {
      const size_t pageIdx = m_ptr / vec_ptr->item_per_page;
      const size_t pageOffset = m_ptr % vec_ptr->item_per_page;
      T* cachePage = (T*)malloc(vec_ptr->item_per_page * sizeof(T));
      if constexpr (AUTH) {
        vec_ptr->server.Read(pageIdx, cachePage, 0);
      } else {
        vec_ptr->server.Read(pageIdx, cachePage);
      }
      T result = *(cachePage + pageOffset);
      free(cachePage);
      return result;
    }

    const T* operator->() { return &(**this); }

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

    page_idx_type get_page_idx() const {
      return m_ptr / vec_ptr->item_per_page;
    }

    page_offset_type get_page_offset() const {
      return m_ptr % vec_ptr->item_per_page;
    }

    pointer get_m_ptr() const { return m_ptr; }

    auto& getVector() { return *vec_ptr; }

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
    pointer m_ptr;
    Vector* vec_ptr;
  };

  // default:
  explicit Vector(
      uint64_t N_, uint64_t item_per_page,
      typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(N_),
        item_per_page(item_per_page),
        server(_backend, N_ / item_per_page, sizeof(T) * item_per_page) {
    Assert(N_ % item_per_page == 0);
  }

  Vector(Vector&& vec)
      : server(vec.server),
        referenceCount(vec.referenceCount),
        N(vec.N),
        item_per_page(vec.item_per_page) {
    vec.server.preventFree();
  }

  Vector(Vector& vec) = delete;

  template <class InputIterator>
  Vector(InputIterator inputBegin, InputIterator inputEnd,
         uint64_t item_per_page,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(inputEnd - inputBegin),
        item_per_page(item_per_page),
        server(_backend, (inputEnd - inputBegin) / item_per_page,
               sizeof(T) * item_per_page) {
    Assert((inputEnd - inputBegin) % item_per_page == 0);
    CopyOut(inputBegin, inputEnd, begin());
  }

  // these functions are slow, as they need to read the whole page
  T& At(uint64_t index) { return *Iterator(index, *this); }
  const T& At(uint64_t index) const { return *Iterator(index, *this); }

  T& operator[](uint64_t index) { return At(index); }
  const T& operator[](uint64_t index) const { return At(index); }

  uint64_t size() const { return N; }

  Iterator begin() { return Iterator(0, *this); }

  Iterator end() { return Iterator(N, *this); }
};

template <class InputIterator, class OutputIterator>
static OutputIterator CopyOut(InputIterator begin, InputIterator end,
                              OutputIterator to, uint32_t counter = 0) {
  using T = typename std::iterator_traits<OutputIterator>::value_type;
  auto& vec = to.getVector();
  using Vec = std::remove_reference_t<decltype(vec)>;
  constexpr bool auth = Vec::auth;
  size_t inputSize = end - begin;
  size_t item_per_page = vec.item_per_page;
  size_t toBeginOffset = to.get_page_offset();
  size_t toBeginPageIdx = to.get_page_idx();
  auto toEnd = to + inputSize;
  Assert(to.getVector().begin() <= to);
  Assert(toEnd <= toEnd.getVector().end());
  size_t toEndOffset = toEnd.get_page_offset();
  size_t toEndPageIdx = toEnd.get_page_idx();

  auto from = begin;
  // printf("item_per_page = %ld, input_size = %ld\n", item_per_page,
  // inputSize);
  Assert(toBeginOffset == 0);
  Assert(toEndOffset == 0);
  Assert(inputSize % item_per_page == 0);
  for (size_t pageIdx = toBeginPageIdx; pageIdx < toEndPageIdx; ++pageIdx) {
    vec.server.Write(pageIdx, &(*from), counter);
    from += item_per_page;
  }

  Assert(from == end);
  return toEnd;
}

template <class InputIterator, class OutputIterator>
static OutputIterator CopyIn(InputIterator begin, InputIterator end,
                             OutputIterator to, uint32_t counter = 0) {
  Assert(end <= end.getVector().end());
  using T = typename std::iterator_traits<OutputIterator>::value_type;
  auto& vec = begin.getVector();
  using Vec = std::remove_reference_t<decltype(vec)>;
  static size_t item_per_page = vec.item_per_page;
  size_t beginOffset = begin.get_page_offset();
  size_t beginPageIdx = begin.get_page_idx();
  size_t endOffset = end.get_page_offset();
  size_t endPageIdx = end.get_page_idx();
  Assert(beginOffset == 0);
  Assert(endOffset == 0);
  for (size_t pageIdx = beginPageIdx; pageIdx < endPageIdx; ++pageIdx) {
    vec.server.Read(pageIdx, &(*to), counter);

    to += item_per_page;
  }

  return to;
}

}  // namespace EM::DynamicPageVector