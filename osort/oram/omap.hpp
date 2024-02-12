#pragma once

#include "oram.hpp"

namespace ODSL {

template <typename K, const short max_fan_out = 9,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct BPlusNode {
  K keys[max_fan_out - 1];
  UidBlock<PositionType, UidType> children[max_fan_out];
  short numChildren = 0;
// friend with cout
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& os, const BPlusNode& node) {
    os << "BPlusNode  ";
    for (int i = 0; i < node.numChildren; i++) {
      os << (int64_t)node.children[i].uid << " @ "
         << (int64_t)node.children[i].data;
      if (i != node.numChildren - 1) {
        os << ", (" << node.keys[i] << "), ";
      }
    }
    return os;
  }
#endif
};

template <typename K, typename V, const short max_chunk_size = 9>
struct BPlusLeaf {
  struct KVPair {
    K key;
    V value;
  };
  KVPair kv[max_chunk_size];
  short numElements = 0;
// friend with cout
#ifndef ENCLAVE_MODE
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
#endif
};

template <typename K, typename V, const short max_fan_out = 9,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct OMap {
  using UidPosition = UidBlock<PositionType, UidType>;
  using BPlusNode_ = BPlusNode<K, max_fan_out, PositionType, UidType>;

  // ensure each leaf contains odd number of elements for performance
  static constexpr short max_chunk_size = std::max(
      3UL, 2 * ((sizeof(BPlusNode_) / (sizeof(K) + sizeof(V))) / 2) + 1);
  using BPlusLeaf_ = BPlusLeaf<K, V, max_chunk_size>;
  using LeafORAM_ = ORAM<BPlusLeaf_, PositionType, UidType>;
  using InternalORAM_ = ORAM<BPlusNode_, PositionType, UidType>;
  LeafORAM_ leafOram;
  std::vector<InternalORAM_> internalOrams;
  std::vector<PositionType> oramSizes;
  PositionType maxSize = 1024;
  static constexpr short min_fan_out = (max_fan_out + 1) / 2;
  static constexpr short min_chunk_size = (max_chunk_size + 1) / 2;

  OMap() {}

  OMap(PositionType maxSize,
       size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) >> 1) {
    SetSize(maxSize, cacheBytes);
  }

  void SetSize(PositionType maxSize,
               size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) >> 1) {
    this->maxSize = maxSize;
    size_t numLevel = 0;
    size_t leafOramSize = divRoundUp(maxSize, min_chunk_size);
    PositionType internalSize = divRoundUp(leafOramSize, min_fan_out);
    for (PositionType size = internalSize;;
         size = divRoundUp(size, min_fan_out)) {
      ++numLevel;
      if (size <= 1) {
        break;
      }
    }
    oramSizes.reserve(numLevel + 1);  // avoid fragmentation
    for (PositionType size = internalSize;;
         size = divRoundUp(size, min_fan_out)) {
      oramSizes.push_back(size);
      if (size <= 1) {
        break;
      }
    }
    internalOrams.reserve(numLevel);
    size_t remainCacheBytes = cacheBytes;
    for (auto it = oramSizes.rbegin(); it != oramSizes.rend(); ++it) {
      PositionType size = *it;
      size_t levelCacheBytes = remainCacheBytes / (numLevel + 1);
      internalOrams.emplace_back(size, levelCacheBytes);
      size_t memUsed = internalOrams.back().GetMemoryUsage();
      remainCacheBytes -= memUsed;
      --numLevel;
    }
    leafOram.SetSize(leafOramSize, remainCacheBytes);
    std::reverse(oramSizes.begin(), oramSizes.end());
    oramSizes.push_back(leafOramSize);
    // printf("OMap: leafOramSize: %lu, internalSize: %lu\n", leafOramSize,
    //        internalSize);
  };

  ~OMap() {}

  template <class Reader>
  void InitFromReader(Reader& reader) {
    printf("InitFromReader\n");
    using T = typename Reader::value_type;
    static_assert(std::is_same_v<T, std::pair<K, V>>,
                  "Reader must read std::pair<K, V>");
    size_t initSize = divRoundUp(reader.size(), max_chunk_size);
    using PositionVec = StdVector<PositionType>;
    using PositionVecWriter = typename PositionVec::Writer;
    using KeyVec = StdVector<K>;
    using KeyWriter = typename KeyVec::Writer;
    PositionVec positions(initSize);
    KeyVec keys(initSize);
    KeyWriter keyWriter(keys.begin(), keys.end());
    /**
      halfFullLeafCount * ((max_chunk_size + 1) / 2) + (initSize -
      halfFullLeafCount) * max_chunk_size >= reader.size(); */
    PositionType halfFullLeafCount =
        (initSize * max_chunk_size - reader.size()) / (max_chunk_size / 2);

    EM::VirtualVector::VirtualReader<BPlusLeaf_> leafReader(
        initSize, [&reader, &keyWriter, halfFullLeafCount](PositionType i) {
          BPlusLeaf_ leaf;
          int j = 0;
          int jUpper =
              i < halfFullLeafCount ? (max_chunk_size + 1) / 2 : max_chunk_size;
          for (; j < jUpper; j++) {
            if (reader.eof()) {
              break;
            }
            auto kv = reader.read();
            leaf.kv[j] = {kv.first, kv.second};
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
    printf("leafOram.InitFromReader done\n");
    // TODO end at a linear oram
    for (auto it = internalOrams.rbegin(); it != internalOrams.rend(); ++it) {
      auto& internalOram = *it;
      PositionType newInitSize = divRoundUp(initSize, max_fan_out);
      PositionType halfFullNodeCount =
          (newInitSize * max_fan_out - initSize) / (max_fan_out / 2);
      PositionVec newPositions(newInitSize);
      KeyVec newKeys(newInitSize);
      KeyWriter newKeyWriter(newKeys.begin(), newKeys.end());
      PositionType idx = 0;
      EM::VirtualVector::VirtualReader<BPlusNode_> nodeReader(
          newInitSize, [&newKeyWriter, &keys, &positions, &idx, initSize,
                        halfFullNodeCount](PositionType i) {
            BPlusNode_ node;
            int j = 0;
            newKeyWriter.write(keys[idx]);
            int jUpper =
                i < halfFullNodeCount ? (max_fan_out + 1) / 2 : max_fan_out;
            for (; j < jUpper; ++j, ++idx) {
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
            Assert(initSize <= max_fan_out || j >= (max_fan_out + 1) / 2);
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
    printf("InitFromReader done\n");
  }

  bool find(const K& key, V& valOut) {
    bool foundFlag = false;
    int level = 0;
    PositionType newPos = 0;
    UidPosition child;
    child.data = 0;
    child.uid = 0;
    for (auto& oram : internalOrams) {
      PositionType childNewPos = UniformRandom(oramSizes[level + 1] - 1);
      UidPosition nextChild;
      std::function<void(BPlusNode_&)> updateFunc = [&](BPlusNode_& node) {
        nextChild = node.children[0];
        short childIdx = 0;
        for (short i = 1; i < max_fan_out; ++i) {
          bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
          obliMove(flag, nextChild, node.children[i]);
          childIdx += flag;
        }
        for (short i = 0; i < max_fan_out; ++i) {
          bool flag = i == childIdx;
          obliMove(flag, node.children[i].data, childNewPos);
        }
      };
      oram.Update(child.data, child.uid, newPos, updateFunc);
      child = nextChild;
      newPos = childNewPos;
      ++level;
    }
    BPlusLeaf_ leaf;
    leafOram.Read(child.data, child.uid, leaf, newPos);
    for (short i = 0; i < max_chunk_size; ++i) {
      bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
      // if (i < leaf.numElements) {
      //   leaf.kv[i].key.print();
      //   printf(" ");
      // }
      foundFlag |= flag;
      obliMove(flag, valOut, leaf.kv[i].value);
    }
    return foundFlag;  // TODO return false if not found
  }

  bool update(const K& key, const std::function<void(V&)>& valUpdateFunc,
              V& valOut) {
    bool foundFlag = false;
    auto leafUpdateFunc = [&](BPlusLeaf_& leaf) {
      short foundIdx = max_chunk_size;
      for (short i = 0; i < max_chunk_size; ++i) {
        bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
        foundFlag |= flag;
        obliMove(flag, valOut, leaf.kv[i].value);
        obliMove(flag, foundIdx, i);
      }
      valUpdateFunc(valOut);
      for (short i = 0; i < max_chunk_size; ++i) {
        bool flag = i == foundIdx;
        obliMove(flag, leaf.kv[i].value, valOut);
      }
    };
    int level = 0;
    PositionType newPos = 0;
    UidPosition child;
    child.data = 0;
    child.uid = 0;
    for (auto& oram : internalOrams) {
      PositionType childNewPos = UniformRandom(oramSizes[level + 1] - 1);
      UidPosition nextChild;
      std::function<void(BPlusNode_&)> updateFunc = [&](BPlusNode_& node) {
        nextChild = node.children[0];
        short childIdx = 0;
        for (short i = 1; i < max_fan_out; ++i) {
          bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
          obliMove(flag, nextChild, node.children[i]);
          childIdx += flag;
        }
        for (short i = 0; i < max_fan_out; ++i) {
          bool flag = i == childIdx;
          obliMove(flag, node.children[i].data, childNewPos);
        }
      };
      oram.Update(child.data, child.uid, newPos, updateFunc);
      child = nextChild;
      newPos = childNewPos;
      ++level;
    }
    leafOram.Update(child.data, child.uid, newPos, leafUpdateFunc);
    return foundFlag;
  }

  bool update(const K& key, const std::function<void(V&)>& valUpdateFunc) {
    V valOut;
    return update(key, valUpdateFunc, valOut);
  }

  bool update(const K& key, const V& val) {
    V valOut;
    auto valUpdateFunc = [&val](V& v) { v = val; };
    return update(key, valUpdateFunc, valOut);
  }

  // returns {splitFlag, foundFlag}
  std::pair<bool, bool> addAndSplitLeaf(BPlusLeaf_& leaf, BPlusLeaf_& newLeaf,
                                        const K& key, const V& val) {
    bool foundFlag = false;
    // update element if exists
    for (short i = 0; i < max_chunk_size; ++i) {
      bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
      foundFlag |= flag;
      obliMove(flag, leaf.kv[i].value, val);
    }
    // add new element if not found and not full
    leaf.numElements += !foundFlag;
    bool splitFlag = leaf.numElements > max_chunk_size;

    // assume split, the index of the rightmost element in the new leaf
    static constexpr short newLeafLastIdx = (max_chunk_size - 1) / 2;
    newLeaf.kv[newLeafLastIdx] = {key, val};

    // insert the new element by swapping from back
    bool lastKVSwapFlag = !splitFlag | (key < leaf.kv[max_chunk_size - 1].key);
    lastKVSwapFlag &= !foundFlag;
    obliSwap(lastKVSwapFlag, leaf.kv[max_chunk_size - 1],
             newLeaf.kv[newLeafLastIdx]);
    for (short i = max_chunk_size - 1; i > 0; --i) {
      bool swapFlag =
          (i >= leaf.numElements) | (leaf.kv[i].key < leaf.kv[i - 1].key);
      swapFlag &= !foundFlag;
      obliSwap(swapFlag, leaf.kv[i - 1], leaf.kv[i]);
    }
    // put the right half of the leaf into the new leaf
    for (short i = max_chunk_size / 2 + 1, j = 0; i < max_chunk_size;
         ++i, ++j) {
      // it's ok to have duplicate kv, since we control validity of
      // element by newElements, and the validity of leaf by uid in oram
      newLeaf.kv[j] = leaf.kv[i];  // TODO optimize with memcpy
    }
    newLeaf.numElements = newLeafLastIdx + 1;
    obliMove(splitFlag, leaf.numElements, (short)(max_chunk_size / 2 + 1));
    // std::cout << "Split flag: " << splitFlag << std::endl;
    // std::cout << "leaf: " << leaf << std::endl;
    // if (splitFlag) std::cout << "newLeaf: " << newLeaf << std::endl;
    return {splitFlag, foundFlag};
  }

  // returns splitFlag
  bool addAndSplitNode(BPlusNode_& node, BPlusNode_& newNode, const K& key,
                       const PositionType& pos, const UidType& uid,
                       short keyRank) {
    Assert(keyRank > 0);  // always insert the right child
    // use keyRank to minimize # of key comparisons

    bool realFlag = uid != DUMMY<UidType>();
    obliMove(!realFlag, keyRank, max_fan_out);  // ensure no swap

    node.numChildren += realFlag;
    bool splitFlag = node.numChildren > max_fan_out;

    // assume split, the index of the rightmost element in the new node
    static constexpr short newNodeLastIdx = (max_fan_out - 1) / 2;
    static_assert(newNodeLastIdx != 0);
    newNode.keys[newNodeLastIdx - 1] = key;
    newNode.children[newNodeLastIdx] = {pos, uid};

    // insert the new element by swapping from back
    bool lastKVSwapFlag = keyRank < max_fan_out;
    obliSwap(lastKVSwapFlag, node.children[max_fan_out - 1],
             newNode.children[newNodeLastIdx]);
    obliSwap(lastKVSwapFlag, node.keys[max_fan_out - 2],
             newNode.keys[newNodeLastIdx - 1]);
    for (short i = max_fan_out - 1; i > 1; --i) {
      bool swapFlag = keyRank < i;
      obliSwap(swapFlag, node.keys[i - 2], node.keys[i - 1]);
      obliSwap(swapFlag, node.children[i - 1], node.children[i]);
    }
    // put the right half of the node into the new node
    for (short i = max_fan_out / 2 + 1, j = 0; i < max_fan_out; ++i, ++j) {
      // it's ok to have duplicate kv, since we control validity of
      // element by newElements, and the validity of node by uid in oram
      newNode.children[j] = node.children[i];  // TODO optimize with memcpy
      if (j > 0) {
        newNode.keys[j - 1] = node.keys[i - 1];
      }
    }

    newNode.numChildren = newNodeLastIdx + 1;
    obliMove(splitFlag, node.numChildren, (short)(max_fan_out / 2 + 1));
    // std::cout << "Split flag: " << splitFlag << std::endl;
    // std::cout << "node: " << node << std::endl;
    // if (splitFlag) std::cout << "newNode: " << newNode << std::endl;
    // node[max_fan_out / 2] will be the key passed to the parent if splitFlag
    return splitFlag;
  }

  // will update rightParentKey if real
  static void redistributeNodeLeftToRight(bool real, BPlusNode_& nodeLeft,
                                          BPlusNode_& nodeRight,
                                          K& rightParentKey) {
    Assert(!real | (nodeRight.numChildren == (max_fan_out - 1) / 2));
    Assert(!real | (nodeLeft.numChildren > (max_fan_out + 1) / 2));
    for (short i = (max_fan_out - 1) / 2 - 1; i > 0; --i) {
      obliMove(real, nodeRight.children[i + 1], nodeRight.children[i]);
      obliMove(real, nodeRight.keys[i], nodeRight.keys[i - 1]);
    }
    obliMove(real, nodeRight.children[1], nodeRight.children[0]);
    obliMove(real, nodeRight.keys[0], rightParentKey);
    for (short i = (max_fan_out + 1) / 2; i < max_fan_out; ++i) {
      bool movFlag = real & (i == nodeLeft.numChildren - 1);
      obliMove(movFlag, nodeRight.children[0], nodeLeft.children[i]);
      obliMove(movFlag, rightParentKey, nodeLeft.keys[i - 1]);
    }
    nodeLeft.numChildren -= (short)real;
    nodeRight.numChildren += (short)real;
  }

  static void redistributeNodeRightToLeft(bool real, BPlusNode_& nodeLeft,
                                          BPlusNode_& nodeRight,
                                          K& rightParentKey) {
    // if right count = left count + 1, we move an entry to left,
    // this makes coalesce easier
    Assert(!real | (nodeRight.numChildren >= (max_fan_out + 1) / 2));
    Assert(!real | (nodeLeft.numChildren == (max_fan_out - 1) / 2));
    Assert(!real | (nodeLeft.numChildren < nodeRight.numChildren));
    short leftNextSlotIdx = (max_fan_out - 1) / 2;
    obliMove(real, nodeLeft.children[leftNextSlotIdx], nodeRight.children[0]);
    obliMove(real, nodeLeft.keys[leftNextSlotIdx - 1], rightParentKey);
    obliMove(real, rightParentKey, nodeRight.keys[0]);
    obliMove(real, nodeRight.children[0], nodeRight.children[1]);
    for (short i = 1; i < max_fan_out - 1; ++i) {
      obliMove(real, nodeRight.children[i], nodeRight.children[i + 1]);
      obliMove(real, nodeRight.keys[i - 1], nodeRight.keys[i]);
    }
    nodeLeft.numChildren += (short)real;
    nodeRight.numChildren -= (short)real;
  }

  // returns if coalesce
  static bool redistributeOrCoalesceNode(bool real, BPlusNode_& nodeLeft,
                                         BPlusNode_& nodeRight,
                                         K& rightParentKey) {
    Assert(!real |
           (nodeLeft.numChildren + nodeRight.numChildren >= max_fan_out));
    bool coalesceFlag =
        real & (nodeLeft.numChildren + nodeRight.numChildren <= max_fan_out);

    bool leftToRightFlag = real &
                           (nodeRight.numChildren < (max_fan_out + 1) / 2) &
                           (nodeLeft.numChildren > (max_fan_out + 1) / 2);
    bool rightToLeftFlag =
        real & (nodeLeft.numChildren < (max_fan_out + 1) / 2);

    // std::cout << "nodeLeft: " << nodeLeft << std::endl;
    // std::cout << "nodeRight: " << nodeRight << std::endl;
    // std::cout << "flags: leftToRightFlag: " << leftToRightFlag
    //           << ", rightToLeftFlag: " << rightToLeftFlag << std::endl
    //           << std::endl;
    redistributeNodeLeftToRight(leftToRightFlag, nodeLeft, nodeRight,
                                rightParentKey);
    redistributeNodeRightToLeft(rightToLeftFlag, nodeLeft, nodeRight,
                                rightParentKey);
    Assert(!coalesceFlag | (nodeLeft.numChildren == (max_fan_out + 1) / 2));
    Assert(!coalesceFlag | (nodeRight.numChildren == max_fan_out / 2));
    for (short i = 0; i < max_fan_out / 2; ++i) {
      obliMove(coalesceFlag, nodeLeft.children[i + (max_fan_out + 1) / 2],
               nodeRight.children[i]);
      if (i > 0) {
        obliMove(coalesceFlag, nodeLeft.keys[i + (max_fan_out + 1) / 2 - 1],
                 nodeRight.keys[i - 1]);
      }
    }
    obliMove(coalesceFlag, nodeLeft.keys[(max_fan_out + 1) / 2 - 1],
             rightParentKey);
    obliMove(coalesceFlag, nodeLeft.numChildren, max_fan_out);
    // std::cout << "After: " << std::endl;
    // std::cout << "nodeLeft: " << nodeLeft << std::endl;
    // std::cout << "nodeRight: " << nodeRight << std::endl;
    // std::cout << "flags: leftToRightFlag: " << leftToRightFlag
    //           << ", rightToLeftFlag: " << rightToLeftFlag << std::endl
    //           << std::endl
    //           << std::endl;
    return coalesceFlag;
  }

  // will update rightParentKey if real
  static void redistributeLeafLeftToRight(bool real, BPlusLeaf_& leafLeft,
                                          BPlusLeaf_& leafRight,
                                          K& rightParentKey) {
    Assert(!real | rightParentKey <= leafRight.kv[0].key);
    Assert(!real | (leafRight.numElements == (max_chunk_size - 1) / 2));
    Assert(!real | (leafLeft.numElements > (max_chunk_size + 1) / 2));
    for (short i = (max_chunk_size - 1) / 2 - 1; i >= 0; --i) {
      obliMove(real, leafRight.kv[i + 1], leafRight.kv[i]);
    }
    for (short i = (max_chunk_size + 1) / 2; i < max_chunk_size; ++i) {
      bool movFlag = real & (i == leafLeft.numElements - 1);
      obliMove(movFlag, leafRight.kv[0], leafLeft.kv[i]);
    }
    rightParentKey = leafRight.kv[0].key;
    leafLeft.numElements -= (short)real;
    leafRight.numElements += (short)real;
  }

  static void redistributeLeafRightToLeft(bool real, BPlusLeaf_& leafLeft,
                                          BPlusLeaf_& leafRight,
                                          K& rightParentKey) {
    // if right count = left count + 1, we move an entry to left,
    // this makes coalesce easier
    Assert(!real | rightParentKey <= leafRight.kv[0].key);
    Assert(!real | (leafRight.numElements >= (max_chunk_size + 1) / 2));
    Assert(!real | (leafLeft.numElements == (max_chunk_size - 1) / 2));
    Assert(!real | (leafLeft.numElements < leafRight.numElements));
    short leftNextSlotIdx = (max_chunk_size - 1) / 2;
    obliMove(real, leafLeft.kv[leftNextSlotIdx], leafRight.kv[0]);
    for (short i = 0; i < max_chunk_size - 1; ++i) {
      obliMove(real, leafRight.kv[i], leafRight.kv[i + 1]);
    }
    rightParentKey = leafRight.kv[0].key;
    leafLeft.numElements += (short)real;
    leafRight.numElements -= (short)real;
  }

  // returns if coalesce
  static bool redistributeOrCoalesceLeaf(bool real, BPlusLeaf_& leafLeft,
                                         BPlusLeaf_& leafRight,
                                         K& rightParentKey) {
    // std::cout << "leafLeft: " << leafLeft << std::endl;
    // std::cout << "leafRight: " << leafRight << std::endl;
    Assert(!real |
           (leafLeft.numElements + leafRight.numElements >= max_chunk_size));
    bool coalesceFlag =
        real & (leafLeft.numElements + leafRight.numElements <= max_chunk_size);

    bool leftToRightFlag = real &
                           (leafRight.numElements < (max_chunk_size + 1) / 2) &
                           (leafLeft.numElements > (max_chunk_size + 1) / 2);
    bool rightToLeftFlag =
        real & (leafLeft.numElements < (max_chunk_size + 1) / 2);
    redistributeLeafLeftToRight(leftToRightFlag, leafLeft, leafRight,
                                rightParentKey);
    redistributeLeafRightToLeft(rightToLeftFlag, leafLeft, leafRight,
                                rightParentKey);
    Assert(!coalesceFlag | (leafLeft.numElements == (max_chunk_size + 1) / 2));
    Assert(!coalesceFlag | (leafRight.numElements == max_chunk_size / 2));
    for (short i = 0; i < max_chunk_size / 2; ++i) {
      obliMove(coalesceFlag, leafLeft.kv[i + (max_chunk_size + 1) / 2],
               leafRight.kv[i]);
    }
    obliMove(coalesceFlag, leafLeft.numElements, max_chunk_size);
    return coalesceFlag;
  }

  // insert a key-value pair into the map, if the key already exist, update the
  // value
  bool insert(const K& key, const V& val) {
    auto oramIt = internalOrams.begin();
    bool foundFlag = false;
    BPlusNode_ newNode;
    PositionType posNewNode;
    bool splitFlag = false;

    K keyNewNode;
    std::function<void(BPlusNode_&)> updateFunc = [&, this](BPlusNode_& node) {
      UidPosition child = node.children[0];
      short childIdx = 0;
      for (short i = 1; i < max_fan_out; ++i) {
        bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);

        // printf("node.kv[%d].key: \n", i);
        // for (int partIdx = 0; partIdx < 5; ++partIdx) {
        //   printf("%X", node.keys[i - 1].part[partIdx]);
        // }
        // printf("\n");

        obliMove(flag, child, node.children[i]);
        childIdx += flag;
      }
      ++oramIt;
      PositionType newPos;
      UidType uidNewNode;
      if (oramIt == internalOrams.end()) {
        BPlusLeaf_ newLeaf;

        // update leaf to write back, and the new node if split
        auto leafUpdateFunc = [&](BPlusLeaf_& leaf) {
          auto sf = addAndSplitLeaf(leaf, newLeaf, key, val);
          splitFlag = sf.first;
          foundFlag = sf.second;
          // for (short i = 0; i < leaf.numElements; ++i) {
          //   printf("leaf.kv[%d].key: \n", i);
          //   for (int partIdx = 0; partIdx < 5; ++partIdx) {
          //     printf("4", leaf.kv[i].key.part[partIdx]);
          //   }
          //   printf("\n");
          // }
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

  // return true if the key is found and erased
  bool eraseEntryFromLeaf(BPlusLeaf_& leaf, const K& key) {
    bool foundFlag = false;
    for (short i = 0; i < max_chunk_size; ++i) {
      bool flag = (i < leaf.numElements) & (key == leaf.kv[i].key);
      foundFlag |= flag;
      if (i != max_chunk_size - 1) {
        obliMove(foundFlag, leaf.kv[i], leaf.kv[i + 1]);
      }
    }
    leaf.numElements -= (short)foundFlag;
    return foundFlag;
  }

  void eraseEntryFromNode(bool real, BPlusNode_& node, short idx) {
    Assert(!real | (idx > 0));
    for (short i = 1; i < max_fan_out - 1; ++i) {
      bool movFlag = real & (i >= idx);
      obliMove(movFlag, node.keys[i - 1], node.keys[i]);
      obliMove(movFlag, node.children[i], node.children[i + 1]);
    }
    node.numChildren -= (short)real;
  }

  bool erase(const K& key) {
    auto oramIt = internalOrams.begin();
    bool foundFlag = false;

    std::function<void(BPlusNode_&)> updateFunc = [&, this](BPlusNode_& node) {
      UidPosition childInfo = node.children[0];

      short childIdx = 0;
      for (short i = 1; i < max_fan_out; ++i) {
        bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);

        // printf("node.kv[%d].key: \n", i);
        // for (int partIdx = 0; partIdx < 5; ++partIdx) {
        //   printf("%X", node.keys[i - 1].part[partIdx]);
        // }
        // printf("\n");

        obliMove(flag, childInfo, node.children[i]);
        childIdx += flag;
      }
      bool noNeighborFlag = node.numChildren == 1;
      bool neighborOnLeftFlag = childIdx == node.numChildren - 1;
      short childNeighborIdx = childIdx + 1;
      obliMove(neighborOnLeftFlag, childNeighborIdx, (short)(childIdx - 1));
      UidPosition childNeighbor = node.children[0];
      for (short i = 1; i < max_fan_out; ++i) {
        obliMove(i == childNeighborIdx, childNeighbor, node.children[i]);
      }
      ++oramIt;
      UidPosition dummyNeighbor;
      dummyNeighbor.data =
          UniformRandom(oramIt != internalOrams.end() ? (oramIt->size() - 1)
                                                      : (leafOram.size() - 1));
      dummyNeighbor.uid = DUMMY<UidType>();  // won't match any real element
      obliMove(noNeighborFlag, childNeighbor, dummyNeighbor);
      short rightChildIdx = childIdx + 1;
      obliMove(neighborOnLeftFlag, rightChildIdx, childIdx);
      // std::cout << "num children = " << node.numChildren << std::endl;
      // std::cout << "childIdx = " << childIdx
      //           << "no neighbor flag = " << noNeighborFlag
      //           << "neighborOnLeftFlag" << neighborOnLeftFlag << std::endl;
      // std::cout << "node: " << node << std::endl;
      PositionType newPos;
      PositionType newPosNeighbor;
      K rightParentKey = node.keys[max_fan_out - 2];
      bool coalesceFlag = false;
      for (short i = 1; i < max_fan_out - 1; ++i) {
        bool flag = i == rightChildIdx;
        obliMove(flag, rightParentKey, node.keys[i - 1]);
      }  // key in parent that separates child and childNeighbor
      if (oramIt == internalOrams.end()) {
        // update leaf to write back, and the new node if split

        auto batchUpdateFunc =
            [&](std::vector<BPlusLeaf_>& leaves) -> std::vector<bool> {
          BPlusLeaf_& child = leaves[0];
          BPlusLeaf_& neighbor = leaves[1];
          foundFlag = eraseEntryFromLeaf(child, key);

          obliSwap(neighborOnLeftFlag, child, neighbor);
          coalesceFlag = redistributeOrCoalesceLeaf(!noNeighborFlag, child,
                                                    neighbor, rightParentKey);
          obliSwap(neighborOnLeftFlag, child, neighbor);
          return {bool(!neighborOnLeftFlag | !coalesceFlag),
                  bool((neighborOnLeftFlag | !coalesceFlag) & !noNeighborFlag)};
        };

        const auto& newPositions = leafOram.BatchUpdate(
            {childInfo.data, childNeighbor.data},
            {childInfo.uid, childNeighbor.uid}, batchUpdateFunc);
        newPos = newPositions[0];
        newPosNeighbor = newPositions[1];
      } else {
        auto batchUpdateFunc =
            [&](std::vector<BPlusNode_>& children) -> std::vector<bool> {
          BPlusNode_& child = children[0];
          updateFunc(child);
          BPlusNode_& neighbor = children[1];
          obliSwap(neighborOnLeftFlag, child, neighbor);
          coalesceFlag = redistributeOrCoalesceNode(!noNeighborFlag, child,
                                                    neighbor, rightParentKey);
          obliSwap(neighborOnLeftFlag, child, neighbor);
          return {bool(!neighborOnLeftFlag | !coalesceFlag),
                  bool((neighborOnLeftFlag | !coalesceFlag) & !noNeighborFlag)};
        };

        const auto& newPositions = oramIt->BatchUpdate(
            {childInfo.data, childNeighbor.data},
            {childInfo.uid, childNeighbor.uid}, batchUpdateFunc);
        newPos = newPositions[0];
        newPosNeighbor = newPositions[1];
      }

      for (short i = 0; i < max_fan_out; ++i) {
        bool flag = i == childIdx;
        obliMove(flag, node.children[i].data, newPos);
        bool flagNeighbor = !noNeighborFlag & (i == childNeighborIdx);
        obliMove(flagNeighbor, node.children[i].data, newPosNeighbor);
      }
      // std::cout << "Right parent key in erase: " << rightParentKey <<
      // std::endl;
      for (short i = 1; i < max_fan_out; ++i) {
        obliMove(!noNeighborFlag & (i == rightChildIdx), node.keys[i - 1],
                 rightParentKey);
      }
      eraseEntryFromNode(coalesceFlag, node, rightChildIdx);
      // redistributeOrCoalesceNode(coalesceFlag, node,
      // node.children[rightChildIdx],
      //                            rightParentKey);

      --oramIt;  // maintain invariant
    };

    oramIt->Update(0, 0, updateFunc);
    return foundFlag;
  }
};
}  // namespace ODSL