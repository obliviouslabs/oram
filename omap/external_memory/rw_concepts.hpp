#include <concepts>
#include <cstddef>

/**
 * @brief This file defines the concepts for reader and writer.
 *
 */

template <class T, class... U>
concept any_of = std::disjunction_v<std::is_same<T, U>...>;

template <typename Reader, typename T = typename Reader::value_type>
concept Readable = requires(Reader reader) {
  { reader.read() } -> any_of<T, const T&>;
  { reader.eof() } -> std::same_as<bool>;
  { reader.size() } -> std::same_as<size_t>;
};

template <typename Writer, typename T = typename Writer::value_type>
concept Writable = requires(Writer writer, T value) {
  { writer.write(value) };
  { writer.eof() } -> std::same_as<bool>;
  { writer.size() } -> std::same_as<size_t>;
  { writer.flush() };
};