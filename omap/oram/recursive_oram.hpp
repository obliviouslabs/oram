#pragma once

#include "oram.hpp"

namespace ODSL {

template <typename T, typename PositionType = uint64_t>
struct RecursiveORAM {
  typedef PositionType UidType;
  static constexpr short fan_out = std::max(64 / (int)sizeof(PositionType), 2);
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

  using InternalORAM = ORAM<InternalNode, PositionType, UidType>;

  std::vector<InternalORAM> internalOrams;

  using LeafORAM = ORAM<LeafNode, PositionType, UidType>;
  LeafORAM leafOram;
  std::vector<PositionType> oramSizes;
  std::vector<UidType> uids;
  std::vector<short> indices;
  PositionType _size;
  Lock _lock;
  bool hasInited = false;

  RecursiveORAM() {}

  RecursiveORAM(PositionType size) { SetSize(size); }

  RecursiveORAM(PositionType size, size_t cacheBytes) {
    SetSize(size, cacheBytes);
  }

  void SetSize(PositionType size,
               size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3UL / 4UL) {
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

  template <typename Reader>
  PositionType InitFromReaderInPlaceHelper(Reader& reader, int level = 0) {
    if (level == oramSizes.size() - 1) {
      LeafNode leafNode;
      for (short i = 0; i < chunk_size; ++i) {
        if (!reader.eof()) {
          leafNode.data[i] = reader.read();
          // std::cout << "leafNode.data[" << i << "] " << leafNode.data[i]
          //           << std::endl;
        } else {
          break;
        }
      }
      UidType uid = leafOram.GetNextUid();
      // printf("leaf write uid %lu\n", uid);
      return leafOram.Write(uid, leafNode);
    }
    InternalNode internalNode;
    for (short i = 0; i < fan_out; ++i) {
      internalNode.children[i] = InitFromReaderInPlaceHelper(reader, level + 1);
      if (reader.eof()) {
        break;
      }
    }
    UidType uid = internalOrams[level].GetNextUid();
    // printf("internal level %d write uid %lu\n", level, uid);
    return internalOrams[level].Write(uid, internalNode);
  }

  template <typename Reader>
  void InitFromReaderInPlace(Reader& reader) {
    if (hasInited) {
      throw std::runtime_error("RecursiveORAM double initialization");
    }
    hasInited = true;
    using ReaderT = typename Reader::value_type;
    static_assert(std::is_same_v<T, ReaderT>, "Reader must reads type T");
    if (reader.size() != _size) {
      throw std::runtime_error("Reader size does not match oram size");
    }

    PositionType rootPos = InitFromReaderInPlaceHelper(reader);
    Assert(rootPos == 0);
  }

  template <class Vec>
  void InitFromVector(Vec& vec) {
    typename Vec::Reader reader(vec.begin(), vec.end());
    InitFromReaderInPlace(reader);
  }

  void InitDefault(const T& defaultValue) {
    EM::VirtualVector::VirtualReader<T> reader(
        _size, [&](PositionType) { return defaultValue; });
    InitFromReaderInPlace(reader);
  }

  template <class Func>
  void Access(UidType address, const Func& accessor) {
    // printf("Access %lu\n", address);
    Critical section(_lock);
    Assert(hasInited);
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
      // printf("level %d: pos %lu, uid %lu, newPos %lu\n", level, pos,
      //        uids[level], newPos);
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
    // printf("level %ld: pos %lu, uid %lu, newPos %lu\n", oramSizes.size() - 1,
    //        pos, uids.back(), newPos);
    leafOram.Update(pos, uids.back(), newPos, updateFunc);
  }

  struct WriteBackBuffer {
    UidType* uids;
    PositionType* newPoses;
    InternalNode* internalNodes;
    LeafNode* leafNodes;
    uint8_t* buffer = NULL;
    size_t bufferSize = 0;
    int numLevel;
    uint64_t batchSize;
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

  WriteBackBuffer oramWriteBackBuffer;

  template <class Func>
  void BatchAccessDeferWriteBack(const std::vector<UidType>& address,
                                 const Func& accessor,
                                 WriteBackBuffer& writeBackBuffer) {
    // printf("Access %lu\n", address);
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

  void WriteBack() { WriteBack(oramWriteBackBuffer); }

  template <class Func>
  void BatchAccessDeferWriteBack(const std::vector<UidType>& address,
                                 const Func& accessor) {
    BatchAccessDeferWriteBack(address, accessor, oramWriteBackBuffer);
  }

  void Read(UidType address, T& out) {
    Access(address, [&](T& data) { out = data; });
  }

  void Write(UidType address, const T& in) {
    Access(address, [&](T& data) { data = in; });
  }

  void BatchReadDeferWriteBack(const std::vector<UidType>& address,
                               std::vector<T>& out,
                               WriteBackBuffer& writeBackBuffer) {
    BatchAccessDeferWriteBack(
        address, [&](const std::vector<T>& data) { out = data; },
        writeBackBuffer);
  }

  void BatchReadDeferWriteBack(const std::vector<UidType>& address,
                               std::vector<T>& out) {
    BatchReadDeferWriteBack(address, out, oramWriteBackBuffer);
  }
};
}  // namespace ODSL