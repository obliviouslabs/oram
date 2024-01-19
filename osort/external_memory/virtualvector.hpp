#pragma once
#include <functional>
namespace EM::VirtualVector {
template <typename VT, class BaseVec>
struct Vector {
  using T = typename BaseVec::value_type;
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  size_t _size;
  std::function<VT(size_t, const T&)> virtualize;
  std::function<T(size_t, const VT&)> devirtualize;
  BaseVec& baseVec;
  Vector(BaseVec& baseVec, std::function<VT(size_t, const T&)> virtualize,
         std::function<T(size_t, const VT&)> devirtualize)
      : baseVec(baseVec),
        _size(baseVec.size()),
        virtualize(virtualize),
        devirtualize(devirtualize) {}

  struct Iterator {
    using iterator_category = std::random_access_iterator_tag;
    using value_type = VT;
    using vector_type = Vector<VT, BaseVec>;
    using difference_type = int64_t;
    using pointer = uint64_t;
    using reference = VT&;
    using const_reference = const VT&;
    constexpr static bool random_access = false;
    Iterator() : m_ptr(0) {}

    // Constructor to initialize the base iterator and custom attribute
    Iterator(pointer ptr, Vector& vec) : m_ptr(ptr), vec_ptr(&vec) {}

    // Copy constructor
    Iterator(const Iterator& it) : m_ptr(it.m_ptr), vec_ptr(it.vec_ptr) {}

    // Copy assignment
    Iterator& operator=(const Iterator& it) {
      m_ptr = it.m_ptr;
      vec_ptr = it.vec_ptr;
      return *this;
    }

    // Destructor
    ~Iterator() {}

    pointer get_m_ptr() const { return m_ptr; }

    Vector* get_vec_ptr() const { return vec_ptr; }

    // Custom functions
    VT operator*() const {
      return vec_ptr->virtualize(m_ptr, vec_ptr->baseVec[m_ptr]);
    }
    pointer operator->() { return m_ptr; }
    const_pointer operator->() const { return m_ptr; }
    Iterator& operator++() {
      ++m_ptr;
      return *this;
    }
    Iterator operator++(int) {
      Iterator tmp(*this);
      ++m_ptr;
      return tmp;
    }
    Iterator& operator--() {
      --m_ptr;
      return *this;
    }
    Iterator operator--(int) {
      Iterator tmp(*this);
      --m_ptr;
      return tmp;
    }
    Iterator& operator+=(int64_t size) {
      m_ptr += size;
      return *this;
    }
    Iterator& operator-=(int64_t size) {
      m_ptr -= size;
      return *this;
    }
    Iterator operator+(int64_t size) const {
      return Iterator(m_ptr + size, *vec_ptr);
    }
    Iterator operator-(int64_t size) const {
      return Iterator(m_ptr - size, *vec_ptr);
    }
    int64_t operator-(const Iterator& it) const { return m_ptr - it.m_ptr; }

   private:
    uint64_t m_ptr;
    Vector* vec_ptr;
  };

  size_t size() const { return _size; }
  Iterator begin() { return Iterator(0, *this); }
  Iterator end() { return Iterator(size(), *this); }

  const_reference operator[](size_t i) const {
    return virtualize(i, baseVec[i]);
  }

  static typename BaseVec::Iterator toBaseIterator(Iterator it) {
    return it.get_vec_ptr()->baseVec.begin() + it.get_m_ptr();
  }

  struct Reader {
    typename BaseVec::Reader baseReader;
    Vector* vec_ptr;
    uint64_t idx = 0;
    Reader() {}
    Reader(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : baseReader(toBaseIterator(_begin), toBaseIterator(_end), _auth),
          vec_ptr(_begin.get_vec_ptr()) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
      baseReader.init(toBaseIterator(_begin), toBaseIterator(_end), _auth);
      vec_ptr = _begin.get_vec_ptr();
    }

    VT get() { return vec_ptr->virtualize(idx, baseReader.get()); }

    VT read() { return vec_ptr->virtualize(idx++, baseReader.read()); }

    bool eof() { return baseReader.eof(); }
  };

  struct PrefetchReader {
    typename BaseVec::PrefetchReader baseReader;
    Vector* vec_ptr;
    uint64_t idx = 0;
    PrefetchReader() {}
    PrefetchReader(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : baseReader(toBaseIterator(_begin), toBaseIterator(_end), _auth),
          vec_ptr(_begin.get_vec_ptr()) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
      baseReader.init(toBaseIterator(_begin), toBaseIterator(_end), _auth);
      vec_ptr = _begin.get_vec_ptr();
    }

    VT get() { return vec_ptr->virtualize(idx, baseReader.get()); }

    VT read() { return vec_ptr->virtualize(idx++, baseReader.read()); }

    bool eof() { return baseReader.eof(); }
  };

  struct Writer {
    typename BaseVec::Writer baseWriter;
    Vector* vec_ptr;
    uint64_t idx = 0;
    Writer() {}
    Writer(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : baseWriter(toBaseIterator(_begin), toBaseIterator(_end), _auth),
          vec_ptr(_begin.get_vec_ptr()) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
      baseWriter.init(toBaseIterator(_begin), toBaseIterator(_end), _auth);
      vec_ptr = _begin.get_vec_ptr();
    }

    void write(const VT& element) {
      return baseWriter.write(vec_ptr->devirtualize(idx++, element));
    }

    bool eof() { return baseWriter.eof(); }

    void flush() { baseWriter.flush(); }
  };
};

template <class InputIterator, class OutputIterator>
static OutputIterator CopyIn(InputIterator begin, InputIterator end,
                             OutputIterator to, uint32_t counter = 0) {
  using Vec = typename InputIterator::vector_type;
  typename Vec::Reader reader(begin, end, counter);
  for (; !reader.eof(); ++to) {
    *to = reader.read();
  }
  return to;
}

template <class InputIterator, class OutputIterator>
static OutputIterator CopyOut(InputIterator begin, InputIterator end,
                              OutputIterator to, uint32_t counter = 0) {
  using Vec = typename OutputIterator::vector_type;
  typename Vec::Writer writer(to, to + (end - begin), counter);
  for (InputIterator from = begin; !writer.eof(); ++from) {
    writer.write(*from);
  }
  writer.eof();
  return to + (end - begin);
}

}  // namespace EM::VirtualVector