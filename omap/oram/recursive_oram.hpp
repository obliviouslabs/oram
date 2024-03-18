#pragma once

#include "adaptive_oram.hpp"

/// @brief This file contains the definition of an oram with recursive position
/// maps.

namespace ODSL {

/**
 * @brief Recursive ORAM stores the position map recursively. Each position map
 * is also an ORAM. The interface of recursive oram is the same as the normal
 * ram.
 *
 * @tparam T The type of the data
 * @tparam PositionType The type of the position, default to uint64_t. If the
 * oram is not very large (< 4e9 elements), however, it is faster to use
 * uint32_t.
 */
template <typename T, typename PositionType = uint64_t>
struct RecursiveORAM {
 private:
  typedef PositionType UidType;
  // Each internal node (i.e. position map node) has fan_out children
  static constexpr short fan_out = std::max(64 / (int)sizeof(PositionType), 2);

  /**
   * @brief Defines the internal node of the position map, stores the positions
   * of its children.
   *
   */
  struct InternalNode {
    PositionType children[fan_out];
#ifndef ENCLAVE_MODE
    friend std::ostream& operator<<(std::ostream& os,
                                    const InternalNode& node) {
      os << "InternalNode{";
      for (int i = 0; i < fan_out; ++i) {
        os << node.children[i] << " ";
      }
      os << "}";
      return os;
    }
#endif
  };

  static constexpr short chunk_size = 1;
  /**
   * @brief Defines the node of the oram that stores the actual data. Each
   * node can hold chunk_size data elements.
   *
   */
  struct LeafNode {
    T data[chunk_size];
#ifndef ENCLAVE_MODE
    friend std::ostream& operator<<(std::ostream& os, const LeafNode& node) {
      os << "LeafNode{";
      for (int i = 0; i < chunk_size; ++i) {
        os << node.data[i] << " ";
      }
      os << "}";
      return os;
    }
#endif
  };

  using InternalORAM = AdaptiveORAM::ORAM<InternalNode, PositionType, UidType>;

  // Stores the position maps. internalOrams[0] is the smallest position map.
  std::vector<InternalORAM> internalOrams;

  using LeafORAM = AdaptiveORAM::ORAM<LeafNode, PositionType, UidType>;
  LeafORAM leafOram;  // Stores the actual data
  // The size of all the orams. oramSizes[0] is the size of the smallest oram.
  // oramSizes.back() is the size of the leaf oram.
  std::vector<PositionType> oramSizes;

  // A global buffer to store the uids at each ORAM level, so that we don't need
  // to allocate memory for each access
  std::vector<UidType> uids;
  // A global buffer to store the indices within the node at each ORAM level.
  std::vector<short> indices;
  // The size of the recursive ORAM
  PositionType _size;
  // The lock to protect the ORAM
  Lock _lock;
  // Whether the ORAM has been initialized
  bool hasInited = false;

  /**
   * @brief Helper function to initialize the ORAM from a reader. Handles a
   * subtree.
   *
   * @tparam Reader The type of the reader
   * @param reader The reader to read the data from
   * @param level The ORAM level of the root of the subtree.
   * @return PositionType The position of the root of the subtree in the oram it
   * belongs to
   */
  template <typename Reader>
  PositionType InitFromReaderHelper(Reader& reader, int level = 0) {
    if (level == oramSizes.size() - 1) {
      LeafNode leafNode;
      for (short i = 0; i < chunk_size; ++i) {
        if (!reader.eof()) {
          leafNode.data[i] = reader.read();
        } else {
          break;
        }
      }
      UidType uid = leafOram.GetNextUid();
      return leafOram.Write(uid, leafNode);
    }
    InternalNode internalNode;
    for (short i = 0; i < fan_out; ++i) {
      internalNode.children[i] = InitFromReaderHelper(reader, level + 1);
      if (reader.eof()) {
        break;
      }
    }
    UidType uid = internalOrams[level].GetNextUid();
    return internalOrams[level].Write(uid, internalNode);
  }

 public:
  RecursiveORAM() {}

  /**
   * @brief Construct a new Recursive ORAM of given size
   *
   * @param size The size of the ORAM
   */
  RecursiveORAM(PositionType size) { SetSize(size); }

  /**
   * @brief Construct a new Recursive ORAM of given size and available cache
   * size
   *
   * @param size The size of the ORAM
   * @param cacheBytes The available cache size in bytes
   */
  RecursiveORAM(PositionType size, size_t cacheBytes) {
    SetSize(size, cacheBytes);
  }

  /**
   * @brief Allocate resource for a Recursive ORAM of given size and available
   * cache size. Note that this function does not fully initialize the ORAM. You
   * need to call InitFromReader or InitDefault to generate the initial position
   * maps and data.
   *
   * @param size The size of the ORAM
   * @param cacheBytes The available cache size in bytes
   */
  void SetSize(PositionType size, size_t cacheBytes = DEFAULT_HEAP_SIZE) {
    _size = size;
    PositionType leafOramSize = divRoundUp(size, chunk_size);
    int numLevel = 0;
    PositionType internalSize = divRoundUp(leafOramSize, fan_out);
    for (PositionType size = internalSize;; size = divRoundUp(size, fan_out)) {
      oramSizes.push_back(size);
      ++numLevel;
      if (size <= 1) {
        break;
      }
    }
    std::reverse(oramSizes.begin(), oramSizes.end());
    internalOrams.reserve(numLevel);
    // Give each level an equal amount of cache at first. Allocate ORAMs from
    // the smallest to the largest. It's likely that the smaller ORAMs do not
    // consume all its budget, so we can tune the budget for the larger ORAMs
    // according to the amount of cache size remained.
    size_t remainCacheBytes = cacheBytes;
    for (PositionType size : oramSizes) {
      size_t levelCacheBytes = remainCacheBytes / (numLevel + 1);

      internalOrams.emplace_back(size, levelCacheBytes);

      size_t memUsed = internalOrams.back().GetMemoryUsage();
      remainCacheBytes -= memUsed;
      --numLevel;
    }
    oramSizes.push_back(leafOramSize);

    leafOram.SetSize(leafOramSize, remainCacheBytes);

    uids.resize(oramSizes.size());
    indices.resize(oramSizes.size());
  }

  /**
   * @brief Get the memory usage of this recursive ORAM
   *
   * @return uint64_t the heap memory usage in bytes
   */
  uint64_t GetMemoryUsage() const {
    uint64_t memUsage = 0;
    for (const auto& oram : internalOrams) {
      memUsage += oram.GetMemoryUsage();
    }
    memUsage += leafOram.GetMemoryUsage();
    return memUsage;
  }

  /**
   * @brief Initialize the ORAM from a reader. The reader must read the same
   * type as the ORAM data.
   *
   * @tparam Reader The type of the reader
   * @param reader The reader to read the data from
   */
  template <typename Reader>
  void InitFromReader(Reader& reader) {
    if (hasInited) {
      throw std::runtime_error("RecursiveORAM double initialization");
    }
    hasInited = true;
    using ReaderT = typename Reader::value_type;
    static_assert(std::is_same_v<T, ReaderT>, "Reader must reads type T");
    if (reader.size() != _size) {
      throw std::runtime_error("Reader size does not match oram size");
    }

    PositionType rootPos = InitFromReaderHelper(reader);
    Assert(rootPos == 0);
  }

  /**
   * @brief Initialize the ORAM with a default value
   *
   * @param defaultValue The default value to initialize the ORAM with
   */
  void InitDefault(const T& defaultValue) {
    EM::VirtualVector::VirtualReader<T> reader(
        _size, [&](PositionType) { return defaultValue; });
    InitFromReader(reader);
  }

  /**
   * @brief Access the ORAM at a given address. The accessor function is called
   * with the data at the address. The accessor function may modify the data in
   * the ORAM address.
   *
   * @tparam Func The type of the accessor function
   * @param address The address to access, must be less than the size of the
   * ORAM
   * @param accessor The accessor function
   */
  template <class Func>
  void Access(UidType address, const Func& accessor) {
    Critical section(_lock);
    Assert(hasInited);
    // first calculate the uid at each ORAM level, and the index within the node
    UidType uid = address / chunk_size;
    short index = address % chunk_size;
    for (int level = oramSizes.size() - 1; level >= 0; --level) {
      uids[level] = uid;
      indices[level] = index;
      index = uid % fan_out;
      uid /= fan_out;
    }
    PositionType pos = 0;
    PositionType newPos = 0;
    for (int level = 0; level < oramSizes.size() - 1; ++level) {
      PositionType nextPos;
      // here we pre-calculate the next position, so that we can update the
      // position map in one go
      PositionType nextNewPos = UniformRandom(oramSizes[level + 1] - 1);
      auto updateFunc = [&](InternalNode& node) -> bool {
        PositionType localNextPos;
        for (short i = 0; i < fan_out; ++i) {
          bool match = i == indices[level];
          obliMove(match, localNextPos, node.children[i]);
          obliMove(match, node.children[i], nextNewPos);
        }
        nextPos = localNextPos;
        return true;
      };
      internalOrams[level].Update(pos, uids[level], newPos, updateFunc);
      pos = nextPos;
      newPos = nextNewPos;
    }
    auto updateFunc = [&](LeafNode& node) -> bool {
      if constexpr (chunk_size == 1) {
        accessor(node.data[0]);
      } else {
        T data;
        short index = indices.back();
        for (short i = 0; i < chunk_size; ++i) {
          obliMove(i == index, data, node.data[i]);
        }
        accessor(data);
        for (short i = 0; i < chunk_size; ++i) {
          obliMove(i == index, node.data[i], data);
        }
      }
      return true;
    };
    // update the actual data
    leafOram.Update(pos, uids.back(), newPos, updateFunc);
  }

  /**
   * @brief A custom data structure to store the data to be written back to the
   * ORAM. This is used to batch the write back operations. To prevent
   * contention in multi-threading environment, we want to reduce the number of
   * mallocs and frees. Therefore, we store all the members in a single buffer.
   *
   */
  struct WriteBackBuffer {
    UidType* uids;                // uids at each level
    PositionType* newPoses;       // new positions at each level
    InternalNode* internalNodes;  // nodes at each internal level
    LeafNode* leafNodes;          // nodes at the leaf level
    uint8_t* buffer = NULL;       // start of the buffer
    size_t bufferSize = 0;        // size of the buffer
    int numLevel;        // number of oram levels, including the leaf level
    uint64_t batchSize;  // batch size

    /**
     * @brief Allocate the buffer and calculate the offset of each member
     *
     * @param numLevel The number of oram levels, including the leaf level
     * @param batchSize The batch size
     */
    void Init(int numLevel, uint64_t batchSize) {
      uint64_t allBatchSize = batchSize * numLevel;
      uint64_t offset = 0;
      uint64_t uidOffset = offset;
      offset += sizeof(UidType) * allBatchSize;
      offset = (offset + 7) / 8 * 8;  // align to 8 bytes
      uint64_t newPosOffset = offset;
      offset += sizeof(PositionType) * allBatchSize;
      offset = (offset + 7) / 8 * 8;  // align to 8 bytes
      uint64_t internalNodeOffset = offset;
      offset += sizeof(InternalNode) * batchSize * (numLevel - 1);
      offset = (offset + 7) / 8 * 8;  // align to 8 bytes
      uint64_t leafNodeOffset = offset;
      offset += sizeof(LeafNode) * batchSize;
      if (bufferSize < offset) {
        if (buffer != NULL) {
          free(buffer);
        }
        buffer = (uint8_t*)malloc(offset);
        bufferSize = offset;
      }
      uids = (UidType*)(buffer + uidOffset);
      newPoses = (PositionType*)(buffer + newPosOffset);
      internalNodes = (InternalNode*)(buffer + internalNodeOffset);
      leafNodes = (LeafNode*)(buffer + leafNodeOffset);
      this->numLevel = numLevel;
      this->batchSize = batchSize;
      memset(newPoses, 0,
             sizeof(PositionType) *
                 batchSize);  // first level of new pos is always 0
    }

    INLINE UidType* GetUids(int level) { return uids + level * batchSize; }

    INLINE PositionType* GetNewPoses(int level) {
      return newPoses + level * batchSize;
    }

    INLINE InternalNode* GetInternalNodes(int level) {
      return internalNodes + level * batchSize;
    }

    INLINE LeafNode* GetLeafNodes() { return leafNodes; }

    ~WriteBackBuffer() {
      if (buffer != NULL) {
        free(buffer);
      }
    }
  };

  // The default write back buffer
  WriteBackBuffer oramWriteBackBuffer;

  /**
   * @brief Perform a batch access to the ORAM. The accessor function is called
   * with a vector of data at the addresses. The accessor function may modify
   * the data in the ORAM addresses. The data is then stored in the write back
   * buffer and the lock is held until WriteBack is called. Addresses must be
   * sorted, and if there are duplicate addresses, the only the updated data of
   * the first duplicate address will be written back.
   *
   * @tparam Func The type of the accessor function
   * @param address The addresses to access, must be sorted.
   * @param accessor The accessor function.
   * @param writeBackBuffer The write back buffer to store the updated data.
   */
  template <class Func>
  void BatchAccessDeferWriteBack(const std::vector<UidType>& address,
                                 const Func& accessor,
                                 WriteBackBuffer& writeBackBuffer) {
    // the overall implementation is similar to single access
    _lock.lock();
    Assert(hasInited);
    int numLevel = oramSizes.size();
    size_t batchSize = address.size();
    writeBackBuffer.Init(numLevel, batchSize);
    std::vector<short> indices(numLevel * batchSize);
    for (size_t i = 0; i < address.size(); ++i) {
      UidType uid = address[i] / chunk_size;
      short index = address[i] % chunk_size;
      for (int level = oramSizes.size() - 1; level >= 0; --level) {
        writeBackBuffer.GetUids(level)[i] = uid;
        indices[level * batchSize + i] = index;
        index = uid % fan_out;
        uid /= fan_out;
      }
    }

    std::vector<PositionType> pos(address.size(), 0);
    std::vector<PositionType> nextPos(address.size());

    for (int level = 0; level < oramSizes.size() - 1; ++level) {
      PositionType* nextNewPos = writeBackBuffer.GetNewPoses(level + 1);
      InternalNode* node = writeBackBuffer.GetInternalNodes(level);
      UidType* uids = writeBackBuffer.GetUids(level);
      for (size_t i = 0; i < address.size(); ++i) {
        nextNewPos[i] = UniformRandom(oramSizes[level + 1] - 1);
      }

      internalOrams[level].BatchReadAndRemove(address.size(), &pos[0], uids,
                                              node);
      uint64_t indexOffset = level * batchSize;
      for (int64_t i = address.size() - 1; i >= 0; --i) {
        // in case of duplicate uid we need to copy what we've updated to
        // maintain consistency scan backward since only the first duplicate
        // uid will be updated

        // we first read the positions of the next level
        for (short j = 0; j < fan_out; ++j) {
          bool match = j == indices[indexOffset + i];
          obliMove(match, nextPos[i], node[i].children[j]);
        }
        // then duplicate what we already updated
        if (i != address.size() - 1) {
          bool copyFlag = uids[i] == uids[i + 1];
          obliMove(copyFlag, node[i], node[i + 1]);
        }
        // and finally update the node to point to the new children
        for (short j = 0; j < fan_out; ++j) {
          bool match = j == indices[indexOffset + i];

          obliMove(match, node[i].children[j], nextNewPos[i]);
        }
      }

      std::swap(pos, nextPos);
    }

    leafOram.BatchReadAndRemove(address.size(), &pos[0],
                                writeBackBuffer.GetUids(numLevel - 1),
                                writeBackBuffer.GetLeafNodes());
    std::vector<T> data(address.size());
    LeafNode* leafNodes = writeBackBuffer.GetLeafNodes();
    uint64_t indexOffset = (numLevel - 1) * batchSize;
    for (size_t i = 0; i < address.size(); ++i) {
      if constexpr (chunk_size == 1) {
        data[i] = leafNodes[i].data[0];
      } else {
        short idx = indices[indexOffset + i];
        for (short j = 0; j < chunk_size; ++j) {
          obliMove(j == idx, data[i], leafNodes[i].data[j]);
        }
      }
    }

    accessor(data);
    for (size_t i = 0; i < address.size(); ++i) {
      if constexpr (chunk_size == 1) {
        leafNodes[i].data[0] = data[i];
      } else {
        short idx = indices[indexOffset + i];
        for (short j = 0; j < chunk_size; ++j) {
          obliMove(j == idx, leafNodes[i].data[j], data[i]);
        }
      }
    }
  }

  /**
   * @brief Write back the data in the write back buffer to the ORAM, and
   * release the lock.
   *
   * @param writeBackBuffer The write back buffer
   */
  void WriteBack(WriteBackBuffer& writeBackBuffer) {
    uint64_t batchSize = writeBackBuffer.batchSize;
    int numLevel = writeBackBuffer.numLevel;
    for (int level = 0; level < numLevel - 1; ++level) {
      internalOrams[level].BatchWriteBack(
          batchSize, writeBackBuffer.GetUids(level),
          writeBackBuffer.GetNewPoses(level),
          writeBackBuffer.GetInternalNodes(level),
          std::vector<bool>(batchSize, true));
    }
    leafOram.BatchWriteBack(batchSize, writeBackBuffer.GetUids(numLevel - 1),
                            writeBackBuffer.GetNewPoses(numLevel - 1),
                            writeBackBuffer.GetLeafNodes(),
                            std::vector<bool>(batchSize, true));
    _lock.unlock();
  }

  /**
   * @brief Write back the data in the global write back buffer to the ORAM, and
   * release the lock.
   *
   */
  void WriteBack() { WriteBack(oramWriteBackBuffer); }

  /**
   * @brief Perform a batch access to the ORAM. The accessor function is called
   * with a vector of data at the addresses. The accessor function may modify
   * the data in the ORAM addresses. The data is then stored in the global write
   * back buffer and the lock is held until WriteBack is called. Addresses must
   * be sorted, and if there are duplicate addresses, the only the updated data
   * of the first duplicate address will be written back.
   *
   * @tparam Func The type of the accessor function
   * @param address The addresses to access, must be sorted.
   * @param accessor The accessor function.
   */
  template <class Func>
  void BatchAccessDeferWriteBack(const std::vector<UidType>& address,
                                 const Func& accessor) {
    BatchAccessDeferWriteBack(address, accessor, oramWriteBackBuffer);
  }

  /**
   * @brief Read the data at a given address
   *
   * @param address The address to read
   * @param out The output data
   */
  void Read(UidType address, T& out) {
    Access(address, [&](T& data) { out = data; });
  }

  /**
   * @brief Write the data at a given address
   *
   * @param address The address to write
   * @param in The input data
   */
  void Write(UidType address, const T& in) {
    Access(address, [&](T& data) { data = in; });
  }

  /**
   * @brief Read a batch of data from the ORAM. The addresses must be sorted.
   * Store the data in the write back buffer. The lock is held until WriteBack
   * is called.
   *
   * @param address The addresses to read from. Must be sorted.
   * @param out The output data
   * @param writeBackBuffer The write back buffer to store the updated data
   */
  void BatchReadDeferWriteBack(const std::vector<UidType>& address,
                               std::vector<T>& out,
                               WriteBackBuffer& writeBackBuffer) {
    BatchAccessDeferWriteBack(
        address, [&](const std::vector<T>& data) { out = data; },
        writeBackBuffer);
  }

  /**
   * @brief Read a batch of data from the ORAM. The addresses must be sorted.
   * Store the data in the global write back buffer. The lock is held until
   * WriteBack is called.
   *
   * @param address The addresses to read from. Must be sorted.
   * @param out The output data
   */
  void BatchReadDeferWriteBack(const std::vector<UidType>& address,
                               std::vector<T>& out) {
    BatchReadDeferWriteBack(address, out, oramWriteBackBuffer);
  }
};
}  // namespace ODSL