#pragma once

#include "oram.hpp"

namespace ORAM {

template <typename K, const ushort max_fan_out = 9,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct BPlusNode {
  K keys[max_fan_out - 1];
  UidBlock<PositionType, UidType> children[max_fan_out];
  ushort numChildren = 0;
  // friend with cout
  friend std::ostream& operator<<(std::ostream& os, const BPlusNode& node) {
    os << "BPlusNode  ";
    for (int i = 0; i < node.numChildren; i++) {
      os << node.children[i].data;
      if (i != node.numChildren - 1) {
        os << ", (" << node.keys[i] << "), ";
      }
    }
    return os;
  }
};

template <typename K, typename V, const ushort max_chunk_size = 9>
struct BPlusLeaf {
  K keys[max_chunk_size];
  V values[max_chunk_size];
  ushort numElements = 0;
  // friend with cout
  friend std::ostream& operator<<(std::ostream& os, const BPlusLeaf& leaf) {
    os << "BPlusLeaf(";
    for (int i = 0; i < leaf.numElements; i++) {
      os << "(" << leaf.keys[i] << ", " << leaf.values[i] << ")";
      if (i != leaf.numElements - 1) {
        os << ", ";
      }
    }
    os << ")";
    return os;
  }
};

template <typename K, typename V, const ushort max_fan_out = 8,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct OMap {
  using UidPosition = UidBlock<PositionType, UidType>;
  using BPlusNode_ = BPlusNode<K, max_fan_out, PositionType, UidType>;

  static constexpr ushort max_chunk_size =
      std::max(2UL, sizeof(BPlusNode_) / (sizeof(K) + sizeof(V)));
  using BPlusLeaf_ = BPlusLeaf<K, V, max_chunk_size>;
  using LeafORAM_ = ORAM<BPlusLeaf_, PositionType, UidType>;
  using InternalORAM_ = ORAM<BPlusNode_, PositionType, UidType>;
  LeafORAM_ leafOram;
  std::vector<InternalORAM_> internalOrams;
  PositionType maxSize = 1024;
  static constexpr ushort min_fan_out = (max_fan_out + 1) / 2;
  static constexpr ushort min_chunk_size = (max_chunk_size + 1) / 2;

  OMap(PositionType maxSize)
      : leafOram(divRoundUp(maxSize, min_chunk_size)), maxSize(maxSize) {
    size_t numLevel = 0;
    PositionType internalSize = divRoundUp(leafOram.size(), min_fan_out);
    for (PositionType size = internalSize;;
         size = divRoundUp(size, min_fan_out)) {
      ++numLevel;
      if (size <= 1) {
        break;
      }
    }
    internalOrams.reserve(numLevel);
    for (PositionType size = internalSize;;
         size = divRoundUp(size, min_fan_out)) {
      internalOrams.emplace_back(size);
      if (size <= 1) {
        break;
      }
    }
  };

  ~OMap() {}

  template <class Reader>
  void InitFromReader(Reader& reader) {
    using T = typename Reader::value_type;
    static_assert(std::is_same_v<T, std::pair<K, V>>,
                  "Reader must read std::pair<K, V>");
    size_t initSize = divRoundUp(reader.size(), max_chunk_size);
    printf("initSize: %lu\n", initSize);
    using PositionVec = StdVector<PositionType>;
    using PositionVecWriter = typename PositionVec::Writer;
    using KeyVec = StdVector<K>;
    using KeyWriter = typename KeyVec::Writer;
    PositionVec positions(initSize);
    KeyVec keys(initSize);
    KeyWriter keyWriter(keys.begin(), keys.end());
    EM::VirtualVector::VirtualReader<BPlusLeaf_> leafReader(
        initSize, [&reader, &keyWriter](uint64_t i) {
          BPlusLeaf_ leaf;
          int j = 0;
          for (; j < max_chunk_size; j++) {
            if (reader.eof()) {
              break;
            }
            auto& [k, v] = reader.read();
            leaf.keys[j] = k;
            leaf.values[j] = v;
          }
          leaf.numElements = j;
          keyWriter.write(leaf.keys[0]);
          // std::cout << "leaf: " << leaf << std::endl;
          return leaf;
        });

    PositionVecWriter posWriter(positions.begin(), positions.end());
    EM::VirtualVector::WrappedWriter<UidPosition, PositionVecWriter>
        wrappedPosWriter(posWriter,
                         [](const UidPosition& uidPos) { return uidPos.data; });
    leafOram.InitFromReader(leafReader, wrappedPosWriter);

    for (InternalORAM_& internalOram : internalOrams) {
      size_t newInitSize = divRoundUp(initSize, max_fan_out);
      PositionVec newPositions(newInitSize);
      KeyVec newKeys(newInitSize);
      KeyWriter newKeyWriter(newKeys.begin(), newKeys.end());
      EM::VirtualVector::VirtualReader<BPlusNode_> nodeReader(
          newInitSize,
          [&newKeyWriter, &keys, &positions, initSize](uint64_t i) {
            BPlusNode_ node;
            int j = 0;
            newKeyWriter.write(keys[i * max_fan_out]);
            for (; j < max_fan_out; j++) {
              size_t idx = i * max_fan_out + j;
              if (idx >= initSize) {
                break;
              }
              if (j != 0) {
                node.keys[j - 1] = keys[idx];
              }
              node.children[j].data = positions[idx];
              node.children[j].uid = idx;
            }
            node.numChildren = j;
            // std::cout << "node: " << node << std::endl;
            return node;
          });
      PositionVecWriter newPosWriter(newPositions.begin(), newPositions.end());
      EM::VirtualVector::WrappedWriter<UidPosition, PositionVecWriter>
          wrappedPosWriter(newPosWriter, [](const UidPosition& uidPos) {
            return uidPos.data;
          });
      internalOram.InitFromReader(nodeReader, wrappedPosWriter);
      initSize = newInitSize;
      std::swap(positions, newPositions);
      std::swap(keys, newKeys);
    }
  }

  bool find(const K& key, V& valOut) {
    auto oramIt = internalOrams.rbegin();
    std::function<void(BPlusNode_&)> updateFunc =
        [&key, &valOut, &oramIt, &updateFunc, this](BPlusNode_& node) {
          UidPosition child = node.children[0];
          short childIdx = 0;
          for (short i = 1; i < max_fan_out; ++i) {
            bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
            obliMove(flag, child, node.children[i]);
            childIdx += flag;
          }
          ++oramIt;
          PositionType newPos;
          // printf("pos: %lu\n", child.data);
          if (oramIt == internalOrams.rend()) {
            BPlusLeaf_ leaf;
            newPos = leafOram.Read(child.data, child.uid, leaf);
            for (short i = 0; i < max_chunk_size; ++i) {
              bool flag = (i < leaf.numElements) & (key == leaf.keys[i]);
              obliMove(flag, valOut, leaf.values[i]);
            }
          } else {
            newPos = oramIt->Update(child.data, child.uid,
                                    updateFunc);  // recursive lambda
          }
          // printf("newPos: %lu\n", newPos);
          for (short i = 0; i < max_fan_out; ++i) {
            bool flag = i == childIdx;
            obliMove(flag, node.children[i].data, newPos);
          }
        };

    oramIt->Update(0, 0, updateFunc);

    return true;  // TODO return false if not found
  }
};
}  // namespace ORAM