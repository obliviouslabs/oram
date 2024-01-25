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
  struct KVPair {
    K key;
    V value;
  };
  KVPair kv[max_chunk_size];
  ushort numElements = 0;
  // friend with cout
  friend std::ostream& operator<<(std::ostream& os, const BPlusLeaf& leaf) {
    os << "BPlusLeaf(";
    for (int i = 0; i < leaf.numElements; i++) {
      os << "(" << leaf.kv[i].key << ", " << leaf.kv[i].value << ")";
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
    // printf("initSize: %lu\n", initSize);
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
            leaf.kv[j] = {k, v};
          }
          leaf.numElements = j;
          keyWriter.write(leaf.kv[0].key);
          // std::cout << "leaf: " << leaf << std::endl;
          return leaf;
        });

    PositionVecWriter posWriter(positions.begin(), positions.end());
    EM::VirtualVector::WrappedWriter<UidPosition, PositionVecWriter>
        wrappedPosWriter(posWriter,
                         [](const UidPosition& uidPos) { return uidPos.data; });
    leafOram.InitFromReader(leafReader, wrappedPosWriter);

    // TODO end at a linear oram
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
    bool foundFlag = false;
    std::function<void(BPlusNode_&)> updateFunc = [&key, &valOut, &oramIt,
                                                   &updateFunc, &foundFlag,
                                                   this](BPlusNode_& node) {
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
          bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
          foundFlag |= flag;
          obliMove(flag, valOut, leaf.kv[i].value);
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

    return foundFlag;  // TODO return false if not found
  }

  // returns {splitFlag, foundFlag}
  std::pair<bool, bool> addAndSplitLeaf(BPlusLeaf_& leaf, BPlusLeaf_& newLeaf,
                                        const K& key, const V& val) {
    bool foundFlag = false;
    // update element if exists
    for (ushort i = 0; i < max_chunk_size; ++i) {
      bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
      foundFlag |= flag;
      obliMove(flag, leaf.kv[i].value, val);
    }
    // add new element if not found and not full
    leaf.numElements += !foundFlag;
    bool splitFlag = leaf.numElements > max_chunk_size;

    // assume split, the index of the rightmost element in the new leaf
    static constexpr ushort newLeafLastIdx = (max_chunk_size - 1) / 2;
    newLeaf.kv[newLeafLastIdx] = {key, val};

    // insert the new element by swapping from back
    bool lastKVSwapFlag = !splitFlag | (key < leaf.kv[max_chunk_size - 1].key);
    lastKVSwapFlag &= !foundFlag;
    obliSwap(lastKVSwapFlag, leaf.kv[max_chunk_size - 1],
             newLeaf.kv[newLeafLastIdx]);
    for (ushort i = max_chunk_size - 1; i > 0; --i) {
      bool swapFlag =
          (i >= leaf.numElements) | (leaf.kv[i].key < leaf.kv[i - 1].key);
      swapFlag &= !foundFlag;
      obliSwap(swapFlag, leaf.kv[i - 1], leaf.kv[i]);
    }
    // put the right half of the leaf into the new leaf
    for (ushort i = max_chunk_size / 2 + 1, j = 0; i < max_chunk_size;
         ++i, ++j) {
      // it's ok to have duplicate kv, since we control validity of
      // element by newElements, and the validity of leaf by uid in oram
      newLeaf.kv[j] = leaf.kv[i];  // TODO optimize with memcpy
    }
    newLeaf.numElements = newLeafLastIdx + 1;
    obliMove(splitFlag, leaf.numElements, (ushort)(max_chunk_size / 2 + 1));
    // std::cout << "Split flag: " << splitFlag << std::endl;
    // std::cout << "leaf: " << leaf << std::endl;
    // if (splitFlag) std::cout << "newLeaf: " << newLeaf << std::endl;
    return {splitFlag, foundFlag};
  }

  // returns splitFlag
  bool addAndSplitNode(BPlusNode_& node, BPlusNode_& newNode, const K& key,
                       const PositionType& pos, const UidType& uid,
                       ushort keyRank) {
    Assert(keyRank > 0);  // always insert the right child
    // use keyRank to minimize # of key comparisons

    bool realFlag = uid != DUMMY<UidType>();
    obliMove(!realFlag, keyRank, max_fan_out);  // ensure no swap

    node.numChildren += realFlag;
    bool splitFlag = node.numChildren > max_fan_out;

    // assume split, the index of the rightmost element in the new node
    static constexpr ushort newNodeLastIdx = (max_fan_out - 1) / 2;
    static_assert(newNodeLastIdx != 0);
    newNode.keys[newNodeLastIdx - 1] = key;
    newNode.children[newNodeLastIdx] = {pos, uid};

    // insert the new element by swapping from back
    bool lastKVSwapFlag = keyRank < max_fan_out;
    obliSwap(lastKVSwapFlag, node.children[max_fan_out - 1],
             newNode.children[newNodeLastIdx]);
    obliSwap(lastKVSwapFlag, node.keys[max_fan_out - 2],
             newNode.keys[newNodeLastIdx - 1]);
    for (ushort i = max_fan_out - 1; i > 1; --i) {
      bool swapFlag = keyRank < i;
      obliSwap(swapFlag, node.keys[i - 2], node.keys[i - 1]);
      obliSwap(swapFlag, node.children[i - 1], node.children[i]);
    }
    // put the right half of the node into the new node
    for (ushort i = max_fan_out / 2 + 1, j = 0; i < max_fan_out; ++i, ++j) {
      // it's ok to have duplicate kv, since we control validity of
      // element by newElements, and the validity of node by uid in oram
      newNode.children[j] = node.children[i];  // TODO optimize with memcpy
      if (j > 0) {
        newNode.keys[j - 1] = node.keys[i - 1];
      }
    }

    newNode.numChildren = newNodeLastIdx + 1;
    obliMove(splitFlag, node.numChildren, (ushort)(max_fan_out / 2 + 1));
    // std::cout << "Split flag: " << splitFlag << std::endl;
    // std::cout << "node: " << node << std::endl;
    // if (splitFlag) std::cout << "newNode: " << newNode << std::endl;
    // node[max_fan_out / 2] will be the key passed to the parent if splitFlag
    return splitFlag;
  }

  // insert a key-value pair into the map, if the key already exist, update the
  // value
  bool insert(const K& key, const V& val) {
    auto oramIt = internalOrams.rbegin();
    bool foundFlag = false;
    BPlusNode_ newNode;
    PositionType posNewNode;
    bool splitFlag = false;

    K keyNewNode;
    std::function<void(BPlusNode_&)> updateFunc = [&, this](BPlusNode_& node) {
      UidPosition child = node.children[0];
      ushort childIdx = 0;
      for (short i = 1; i < max_fan_out; ++i) {
        bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
        obliMove(flag, child, node.children[i]);
        childIdx += flag;
      }
      ++oramIt;
      PositionType newPos;
      UidType uidNewNode;
      if (oramIt == internalOrams.rend()) {
        BPlusLeaf_ newLeaf;

        // update leaf to write back, and the new node if split
        auto leafUpdateFunc = [&](BPlusLeaf_& leaf) {
          auto sf = addAndSplitLeaf(leaf, newLeaf, key, val);
          splitFlag = sf.first;
          foundFlag = sf.second;
        };

        newPos = leafOram.Update(child.data, child.uid, leafUpdateFunc);

        // DUMMY if splitFlag is false
        uidNewNode = leafOram.GetNextUid(splitFlag);
        posNewNode = leafOram.Write(uidNewNode, newLeaf);
        keyNewNode = newLeaf.kv[0].key;
      } else {
        newPos = oramIt->Update(child.data, child.uid,
                                updateFunc);  // recursive lambda
        uidNewNode = oramIt->GetNextUid(splitFlag);
        posNewNode = oramIt->Write(uidNewNode, newNode);
      }

      // printf("newPos: %lu\n", newPos);
      for (short i = 0; i < max_fan_out; ++i) {
        bool flag = i == childIdx;
        obliMove(flag, node.children[i].data, newPos);
      }

      // insert the new node if split
      splitFlag = addAndSplitNode(node, newNode, keyNewNode, posNewNode,
                                  uidNewNode, childIdx + 1);
      keyNewNode = node.keys[max_fan_out / 2];

      --oramIt;  // maintain invariant
    };

    oramIt->Update(0, 0, updateFunc);
    Assert(!splitFlag);  // root should not split
    return foundFlag;
  }
};
}  // namespace ORAM