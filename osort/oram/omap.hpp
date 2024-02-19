#pragma once

#include "oram.hpp"
/**
 * @brief An oblivious map based on B+ tree.
 * Each level of the tree is stored in a separate ORAM.
 * The number of levels is fixed. If a level has only one node, then the node
 * can contain any number of entries, otherwise, the node must contain at least
 * min_fan_out entries.
 */
namespace ODSL {

template <typename K, const short max_fan_out = 9,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct BPlusNode {
  K keys[max_fan_out - 1];  // does not store the key of the first child
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

/**
 * @tparam K must overload < and == operators, should be trivially copyable
 * @tparam V should be trivially copyable
 * @tparam max_fan_out should be odd for best performance
 * @tparam PositionType
 * @tparam UidType
 */
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
       size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3UL / 4UL) {
    SetSize(maxSize, cacheBytes);
  }

  void SetSize(PositionType maxSize,
               size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3UL / 4UL) {
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
    // Each oram gets equal share of the remaining memory. Start from the
    // smallest oram. Deduct the memory used by the oram from the remaining
    // memory.
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

  void Init() {
    BPlusLeaf_ emptyLeaf;
    emptyLeaf.numElements = 0;
    UidType uid = leafOram.GetNextUid();
    PositionType pos = leafOram.Write(uid, emptyLeaf);
    for (auto oramIt = internalOrams.rbegin(); oramIt != internalOrams.rend();
         ++oramIt) {
      BPlusNode_ node;
      node.numChildren = 1;
      node.children[0].data = pos;
      node.children[0].uid = uid;
      uid = oramIt->GetNextUid();
      pos = oramIt->Write(uid, node);
    }
    Assert(uid == 0);
  }

  /**
   * @brief Initialize the oram, the input must be sorted by key.
   * Build the orams bottom up.
   *
   * @tparam Reader
   * @param reader reads std::pair<K, V>
   */
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
    /** satisfies the following condition:
      halfFullLeafCount * (min_chunk_size) + (initSize -
      halfFullLeafCount) * max_chunk_size >= reader.size(); */
    PositionType halfFullLeafCount =
        (initSize * max_chunk_size - reader.size()) / (max_chunk_size / 2);

    EM::VirtualVector::VirtualReader<BPlusLeaf_> leafReader(
        initSize, [&reader, &keyWriter, halfFullLeafCount](PositionType i) {
          BPlusLeaf_ leaf;
          int j = 0;
          int jUpper = i < halfFullLeafCount ? min_chunk_size : max_chunk_size;
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
            int jUpper = i < halfFullNodeCount ? min_fan_out : max_fan_out;
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
            Assert(initSize <= max_fan_out || j >= min_fan_out);
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

  template <class Reader>
  void InitFromReaderInPlace(Reader& reader) {
    using T = typename Reader::value_type;
    static_assert(std::is_same_v<T, std::pair<K, V>>,
                  "Reader must read std::pair<K, V>");
    PositionType initSize = divRoundUp(reader.size(), max_chunk_size);
    // printf("max chunk size = %d\n", max_chunk_size);
    // printf("init size = %lu\n", initSize);
    PositionType halfFullLeafCount =
        (initSize * max_chunk_size - reader.size()) / (max_chunk_size / 2);
    // printf("half full leaf count = %lu\n", halfFullLeafCount);
    std::vector<BPlusNode_> internalNodeCache(internalOrams.size());
    std::vector<K> internalBeginKeys(internalOrams.size());
    std::vector<PositionType> halfFullNodeCounts(internalOrams.size());
    PositionType levelInitSize = initSize;
    for (int i = internalOrams.size() - 1; i >= 0; --i) {
      PositionType upperLevelInitSize = divRoundUp(levelInitSize, max_fan_out);
      halfFullNodeCounts[i] =
          (upperLevelInitSize * max_fan_out - levelInitSize) /
          (max_fan_out / 2);
      levelInitSize = upperLevelInitSize;
    }
    for (PositionType i = 0; i < initSize; ++i) {
      int jUpper = i < halfFullLeafCount ? min_chunk_size : max_chunk_size;
      int j = 0;
      BPlusLeaf_ leaf;
      for (; j < jUpper; ++j) {
        if (reader.eof()) {
          break;
        }
        auto kv = reader.read();
        leaf.kv[j] = {kv.first, kv.second};
      }
      leaf.numElements = j;
      bool endFlag = i == initSize - 1;
      Assert(endFlag || !reader.eof());
      UidType childUid = leafOram.GetNextUid();
      PositionType childPos = leafOram.Write(childUid, leaf);
      // std::cout << "write leaf: " << leaf << " at pos " << childPos
      //           << " with uid " << childUid << std::endl;
      K childKey = leaf.kv[0].key;
      K nextLevelChildKey;
      for (int j = internalOrams.size() - 1; j >= 0; --j) {
        BPlusNode_& node = internalNodeCache[j];
        int kUpper = max_fan_out;
        if (halfFullNodeCounts[j]) {
          kUpper = min_fan_out;
          --halfFullNodeCounts[j];
        }
        node.children[node.numChildren] = {childPos, childUid};
        if (node.numChildren) {
          node.keys[node.numChildren - 1] = childKey;
        } else {
          internalBeginKeys[j] = childKey;
        }
        ++node.numChildren;
        if (node.numChildren == kUpper || endFlag) {
          childUid = internalOrams[j].GetNextUid();
          childPos = internalOrams[j].Write(childUid, node);
          // std::cout << "write node: " << node << " to level " << j << " at
          // pos "
          //           << childPos << std::endl;
          node.numChildren = 0;
          childKey = internalBeginKeys[j];
        } else {
          break;
        }
      }
    }
    // PositionType childPos
  }

  /**
   * @brief Find by key.
   *
   * @param key
   * @param valOut
   * @return true if key found
   * @return false if key not found
   */
  bool find(const K& key, V& valOut, bool isDummy = false) {
    bool foundFlag = false;
    int level = 0;
    PositionType newPos = 0;
    UidPosition child;
    child.data = 0;
    child.uid = 0;
    auto prevOramIt = internalOrams.begin();
    for (auto oramIt = internalOrams.begin(); oramIt != internalOrams.end();
         ++oramIt) {
      auto& oram = *oramIt;
      PositionType childORAMSize = oramSizes[level + 1];
      PositionType childNewPos = UniformRandom(childORAMSize - 1);
      UidPosition nextChild;
      std::function<bool(BPlusNode_&)> updateFunc =
          [&](BPlusNode_& node) -> bool {
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
        return true;
      };
      obliMove(isDummy, child.data,
               UniformRandom(childORAMSize - 1));  // set a random read path
      obliMove(isDummy, child.uid, DUMMY<UidType>());
      oram.lock();
      if (level > 0) {
        prevOramIt->unlock();
      }
      prevOramIt = oramIt;
      oram.Update(child.data, child.uid, newPos, updateFunc);
      child = nextChild;
      newPos = childNewPos;
      ++level;
    }
    BPlusLeaf_ leaf;
    uint64_t leafPosRand;
    obliMove(isDummy, child.data, UniformRandom(leafOram.size() - 1));
    obliMove(isDummy, child.uid, DUMMY<UidType>());
    leafOram.lock();
    prevOramIt->unlock();
    leafOram.Read(child.data, child.uid, leaf, newPos);
    leafOram.unlock();
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

  /**
   * @brief Update the value by the key.
   *
   * @param key
   * @param valUpdateFunc updates the value in place
   * @param valOut allows reading the updated value
   * @return true if key found
   * @return false if key not found
   */
  bool update(const K& key, const std::function<void(V&)>& valUpdateFunc,
              V& valOut, bool isDummy = false) {
    bool foundFlag = false;
    std::function<bool(BPlusLeaf_&)> leafUpdateFunc =
        [&](BPlusLeaf_& leaf) -> bool {
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
      return true;
    };
    int level = 0;
    PositionType newPos = 0;
    UidPosition child;
    child.data = 0;
    child.uid = 0;
    auto prevOramIt = internalOrams.begin();
    for (auto oramIt = internalOrams.begin(); oramIt != internalOrams.end();
         ++oramIt) {
      auto& oram = *oramIt;
      PositionType childORAMSize = oramSizes[level + 1];
      PositionType childNewPos = UniformRandom(childORAMSize - 1);
      UidPosition nextChild;
      std::function<bool(BPlusNode_&)> updateFunc =
          [&](BPlusNode_& node) -> bool {
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
        return true;
      };

      obliMove(isDummy, child.data,
               UniformRandom(childORAMSize - 1));  // set a random read path
      obliMove(isDummy, child.uid, DUMMY<UidType>());
      oram.lock();
      if (level > 0) {
        prevOramIt->unlock();
      }
      prevOramIt = oramIt;
      oram.Update(child.data, child.uid, newPos, updateFunc);
      child = nextChild;
      newPos = childNewPos;
      ++level;
    }

    obliMove(isDummy, child.data, UniformRandom(leafOram.size() - 1));
    obliMove(isDummy, child.uid, DUMMY<UidType>());
    leafOram.lock();
    prevOramIt->unlock();
    leafOram.Update(child.data, child.uid, newPos, leafUpdateFunc);
    leafOram.unlock();
    return foundFlag;
  }

  bool update(const K& key, const std::function<void(V&)>& valUpdateFunc,
              bool isDummy = false) {
    V valOut;
    return update(key, valUpdateFunc, valOut, isDummy);
  }

  bool update(const K& key, const V& val, bool isDummy = false) {
    V valOut;
    auto valUpdateFunc = [&val](V& v) { v = val; };
    return update(key, valUpdateFunc, valOut, isDummy);
  }

  /**
   * @brief helper function of insert
   *
   * @param leaf leaf that may be split
   * @param newLeaf the new leaf if split, can be anything if not split
   * @param key the key to insert
   * @param val the value to insert
   * @return std::pair<bool, bool> {splitFlag, foundFlag}
   */
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

  /**
   * @brief Helper function for erase. Conditionally put the rightmost entry of
   * nodeLeft to node right, and get the new key in the parent node.
   *
   * @param real if the operation is not dummy
   * @param nodeLeft node on the left (has smaller keys)
   * @param nodeRight node on the right (has larger keys)
   * @param rightParentKey the key between nodeLeft and nodeRight stored in
   * parent, may be updated
   */
  static void redistributeNodeLeftToRight(bool real, BPlusNode_& nodeLeft,
                                          BPlusNode_& nodeRight,
                                          K& rightParentKey) {
    Assert(!real | (nodeRight.numChildren == (max_fan_out - 1) / 2));
    Assert(!real | (nodeLeft.numChildren > min_fan_out));
    for (short i = (max_fan_out - 1) / 2 - 1; i > 0; --i) {
      obliMove(real, nodeRight.children[i + 1], nodeRight.children[i]);
      obliMove(real, nodeRight.keys[i], nodeRight.keys[i - 1]);
    }
    obliMove(real, nodeRight.children[1], nodeRight.children[0]);
    obliMove(real, nodeRight.keys[0], rightParentKey);
    for (short i = min_fan_out; i < max_fan_out; ++i) {
      bool movFlag = real & (i == nodeLeft.numChildren - 1);
      obliMove(movFlag, nodeRight.children[0], nodeLeft.children[i]);
      obliMove(movFlag, rightParentKey, nodeLeft.keys[i - 1]);
    }
    nodeLeft.numChildren -= (short)real;
    nodeRight.numChildren += (short)real;
  }

  /**
   * @brief Similar to redistributeNodeLeftToRight, but in the opposite
   * direction
   */
  static void redistributeNodeRightToLeft(bool real, BPlusNode_& nodeLeft,
                                          BPlusNode_& nodeRight,
                                          K& rightParentKey) {
    // if right count = left count + 1, we move an entry to left,
    // this makes coalesce easier
    Assert(!real | (nodeRight.numChildren >= min_fan_out));
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

  /**
   * @brief Redistribute or coalesce the nodes depending on the loads of the
   * nodes. If coalesce, put all entries in the right node into the left node.
   *
   * @param real if the operation is not dummy
   * @param nodeLeft
   * @param nodeRight
   * @param rightParentKey
   * @return true if coalesce
   * @return false if not coalesce
   */
  static bool redistributeOrCoalesceNode(bool real, BPlusNode_& nodeLeft,
                                         BPlusNode_& nodeRight,
                                         K& rightParentKey) {
    Assert(!real |
           (nodeLeft.numChildren + nodeRight.numChildren >= max_fan_out));
    bool coalesceFlag =
        real & (nodeLeft.numChildren + nodeRight.numChildren <= max_fan_out);

    bool leftToRightFlag = real & (nodeRight.numChildren < min_fan_out) &
                           (nodeLeft.numChildren > min_fan_out);
    bool rightToLeftFlag = real & (nodeLeft.numChildren < min_fan_out);

    // std::cout << "nodeLeft: " << nodeLeft << std::endl;
    // std::cout << "nodeRight: " << nodeRight << std::endl;
    // std::cout << "flags: leftToRightFlag: " << leftToRightFlag
    //           << ", rightToLeftFlag: " << rightToLeftFlag << std::endl
    //           << std::endl;
    redistributeNodeLeftToRight(leftToRightFlag, nodeLeft, nodeRight,
                                rightParentKey);
    redistributeNodeRightToLeft(rightToLeftFlag, nodeLeft, nodeRight,
                                rightParentKey);
    Assert(!coalesceFlag | (nodeLeft.numChildren == min_fan_out));
    Assert(!coalesceFlag | (nodeRight.numChildren == max_fan_out / 2));
    for (short i = 0; i < max_fan_out / 2; ++i) {
      obliMove(coalesceFlag, nodeLeft.children[i + min_fan_out],
               nodeRight.children[i]);
      if (i > 0) {
        obliMove(coalesceFlag, nodeLeft.keys[i + min_fan_out - 1],
                 nodeRight.keys[i - 1]);
      }
    }
    obliMove(coalesceFlag, nodeLeft.keys[min_fan_out - 1], rightParentKey);
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

  // similar operations for leaves
  static void redistributeLeafLeftToRight(bool real, BPlusLeaf_& leafLeft,
                                          BPlusLeaf_& leafRight,
                                          K& rightParentKey) {
    Assert(!real | rightParentKey <= leafRight.kv[0].key);
    Assert(!real | (leafRight.numElements == (max_chunk_size - 1) / 2));
    Assert(!real | (leafLeft.numElements > min_chunk_size));
    for (short i = (max_chunk_size - 1) / 2 - 1; i >= 0; --i) {
      obliMove(real, leafRight.kv[i + 1], leafRight.kv[i]);
    }
    for (short i = min_chunk_size; i < max_chunk_size; ++i) {
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
    Assert(!real | (leafRight.numElements >= min_chunk_size));
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

    bool leftToRightFlag = real & (leafRight.numElements < min_chunk_size) &
                           (leafLeft.numElements > min_chunk_size);
    bool rightToLeftFlag = real & (leafLeft.numElements < min_chunk_size);
    redistributeLeafLeftToRight(leftToRightFlag, leafLeft, leafRight,
                                rightParentKey);
    redistributeLeafRightToLeft(rightToLeftFlag, leafLeft, leafRight,
                                rightParentKey);
    Assert(!coalesceFlag | (leafLeft.numElements == min_chunk_size));
    Assert(!coalesceFlag | (leafRight.numElements == max_chunk_size / 2));
    for (short i = 0; i < max_chunk_size / 2; ++i) {
      obliMove(coalesceFlag, leafLeft.kv[i + min_chunk_size], leafRight.kv[i]);
    }
    obliMove(coalesceFlag, leafLeft.numElements, max_chunk_size);
    return coalesceFlag;
  }

  /**
   * @brief insert a key-value pair into the map, if the key already exist,
   update the value
   *
   * @param key
   * @param val
   * @return true if key found
   * @return false if key not found
   */
  bool insert(const K& key, const V& val, bool isDummy = false) {
    auto oramIt = internalOrams.begin();
    bool foundFlag = false;
    BPlusNode_ newNode;
    PositionType posNewNode;
    bool splitFlag = false;

    K keyNewNode;
    std::function<bool(BPlusNode_&)> updateFunc =
        [&, this](BPlusNode_& node) -> bool {
      UidPosition child = node.children[0];
      short childIdx = 0;
      for (short i = 1; i < max_fan_out; ++i) {
        bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
        obliMove(flag, child, node.children[i]);
        childIdx += flag;
      }
      ++oramIt;
      PositionType newPos;
      UidType uidNewNode;
      obliMove(isDummy, child.uid, DUMMY<UidType>());
      if (oramIt == internalOrams.end()) {
        BPlusLeaf_ newLeaf;

        // update leaf to write back, and the new node if split
        std::function<bool(BPlusLeaf_&)> leafUpdateFunc =
            [&](BPlusLeaf_& leaf) -> bool {
          auto sf = addAndSplitLeaf(leaf, newLeaf, key, val);
          splitFlag = sf.first;
          foundFlag = sf.second;
          return true;
        };
        obliMove(isDummy, child.data,
                 UniformRandom(leafOram.size() - 1));  // set a random read path
        newPos = leafOram.Update(child.data, child.uid, leafUpdateFunc);

        // DUMMY if splitFlag is false
        uidNewNode = leafOram.GetNextUid(splitFlag & !isDummy);
        posNewNode = leafOram.Write(uidNewNode, newLeaf);
        keyNewNode = newLeaf.kv[0].key;
      } else {
        obliMove(isDummy, child.data,
                 UniformRandom(oramIt->size() - 1));  // set a random read path
        newPos = oramIt->Update(child.data, child.uid,
                                updateFunc);  // recursive lambda
        uidNewNode = oramIt->GetNextUid(splitFlag & !isDummy);
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
      return true;
    };
    UidType uid = 0;
    obliMove(isDummy, uid, DUMMY<UidType>());
    oramIt->Update(0, uid, updateFunc);
    Assert(!splitFlag);  // root should not split
    return foundFlag;
  }

  /**
   * @brief erase the key from a leaf
   *
   * @param leaf
   * @param key
   * @return true if the key is found and erased
   * @return false otherwise
   */
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

  /**
   * @brief Erase an entry from an internal node by child idx
   *
   * @param real if the operation is not dummy
   * @param node the node to erase from
   * @param idx > 0 if real (because we only delete children on the right)
   */
  void eraseEntryFromNode(bool real, BPlusNode_& node, short idx) {
    Assert(!real | (idx > 0));
    for (short i = 1; i < max_fan_out - 1; ++i) {
      bool movFlag = real & (i >= idx);
      obliMove(movFlag, node.keys[i - 1], node.keys[i]);
      obliMove(movFlag, node.children[i], node.children[i + 1]);
    }
    node.numChildren -= (short)real;
  }

  /**
   * @brief Erase the key from the map
   *
   * @param key
   * @return true if key found and erased
   * @return false otherwise
   */
  bool erase(const K& key, bool isDummy = false) {
    auto oramIt = internalOrams.begin();
    bool foundFlag = false;

    std::function<bool(BPlusNode_&)> updateFunc =
        [&, this](BPlusNode_& node) -> bool {
      UidPosition childInfo = node.children[0];

      short childIdx = 0;
      for (short i = 1; i < max_fan_out; ++i) {
        bool flag = (i < node.numChildren) & !(key < node.keys[i - 1]);
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
        obliMove(isDummy, childInfo.data, UniformRandom(leafOram.size() - 1));
        obliMove(isDummy, childInfo.uid, DUMMY<UidType>());
        obliMove(isDummy, childNeighbor.data,
                 UniformRandom(leafOram.size() - 1));
        obliMove(isDummy, childNeighbor.uid, DUMMY<UidType>());

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
        obliMove(isDummy, childInfo.data,
                 UniformRandom(oramIt->size() - 1));  // set a random read path
        obliMove(isDummy, childInfo.uid, DUMMY<UidType>());
        obliMove(isDummy, childNeighbor.data,
                 UniformRandom(oramIt->size() - 1));
        obliMove(isDummy, childNeighbor.uid, DUMMY<UidType>());

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
      --oramIt;  // maintain invariant
      return true;
    };
    UidType uid = 0;
    obliMove(isDummy, uid, DUMMY<UidType>());
    oramIt->Update(0, uid, updateFunc);
    return foundFlag;
  }
};
}  // namespace ODSL