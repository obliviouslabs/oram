#pragma once
#include <functional>
/// @brief Virtual vector that allows to access a virtualized memory
namespace EM::VirtualVector {
template <typename VT>
struct Vector {
 private:
  using value_type = VT;
  using reference = VT&;
  using const_reference = const VT&;
  using pointer = VT*;
  using const_pointer = const VT*;
  size_t _size;
  std::function<VT&(size_t)> virtualize;

 public:
  Vector(size_t size, std::function<VT&(size_t)> virtualize)
      : _size(size), virtualize(virtualize) {}
  struct Iterator {
    using iterator_category = std::random_access_iterator_tag;
    using value_type = VT;
    using vector_type = Vector<VT>;
    using difference_type = int64_t;
    using pointer = uint64_t;
    using reference = VT&;
    using const_reference = const VT&;
    Iterator() : m_ptr(0), vec_ptr(NULL) {}

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
    VT& operator*() const { return vec_ptr->virtualize(m_ptr); }
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

  const_reference operator[](size_t i) const { return virtualize(i); }

  struct Reader {
    using value_type = VT;
    Iterator it;
    Iterator end;
    Reader() {}
    Reader(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
      it = _begin;
      end = _end;
    }

    VT get() { return *it; }

    VT read() {
      VT res = *it;
      ++it;
      return res;
    }

    bool eof() { return it.get_m_ptr() == end.get_m_ptr(); }

    size_t size() { return end - it; }
  };

  struct Writer {
    using value_type = VT;
    Iterator it;
    Iterator end;
    Writer() {}
    Writer(Iterator _begin, Iterator _end, uint32_t _auth = 0)
        : it(_begin), end(_end) {}

    void init(Iterator _begin, Iterator _end, uint32_t _auth = 0) {
      it = _begin;
      end = _end;
    }

    void write(const VT& element) {
      *it = element;
      ++it;
    }

    bool eof() { return it.get_m_ptr() == end.get_m_ptr(); }

    void flush() {}

    size_t size() { return end - it; }
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

/// @brief A reader based on a customized virtualization function.
/// @tparam VT The value type the reader reads
template <typename VT>
struct VirtualReader {
  using value_type = VT;
  using iterator_type = uint64_t;
  iterator_type it = 0;
  size_t _size;
  std::function<VT(iterator_type)> virtualize;
  VirtualReader(size_t size, std::function<VT(iterator_type)> virtualize)
      : _size(size), virtualize(virtualize) {}

  VT get() { return virtualize(it); }

  VT read() { return virtualize(it++); }

  bool eof() { return it >= _size; }

  size_t size() { return _size; }

  iterator_type getIterator() { return it; }
};

/// @brief A writer that calls a customized function when performing write.
/// @tparam VT The value type the writer writes
template <typename VT>
struct VirtualWriter {
  using value_type = VT;
  using iterator_type = uint64_t;
  iterator_type it = 0;
  size_t _size;
  std::function<void(iterator_type, const VT&)> devirtualize;
  VirtualWriter(size_t size,
                std::function<void(iterator_type, const VT&)> devirtualize)
      : _size(size), devirtualize(devirtualize) {}

  void write(const VT& element) { devirtualize(it++, element); }

  bool eof() { return it >= _size; }

  size_t size() { return _size; }

  void flush() {}
};

}  // namespace EM::VirtualVector