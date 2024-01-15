#pragma once
#include "bucket.hpp"

namespace ORAM::PathORAM {

template <typename T, const int Z = 4, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct Node {
  Bucket<T, Z, PositionType, UidType> bucket;
  Node* left = nullptr;
  Node* right = nullptr;
  Node() = default;
  Node(const Node& other) = default;
  Node(Node&& other) = default;
  Node& operator=(const Node& other) = default;
  Node& operator=(Node&& other) = default;
};

};  // namespace ORAM::PathORAM
