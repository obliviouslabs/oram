#pragma once

#include "common/utils.hpp"
#include "external_memory/server/serverFrontend.hpp"
#include "external_memory/virtualvector.hpp"
#include "oram_common.hpp"

namespace ODSL {

template <typename T, typename UidType = uint64_t,
          const uint64_t page_size = 4096>
struct Page {
  static constexpr uint64_t item_per_page =
      divRoundUp(page_size, sizeof(UidType) + sizeof(T));
  UidType uids[item_per_page];
  T data[item_per_page];
  using Encrypted_t = FreshEncrypted<Page>;
  Page() {
    for (UidType i = 0; i < item_per_page; i++) {
      uids[i] = DUMMY<UidType>();
    }
  }
};

template <typename T, typename UidType = uint64_t,
          typename PageIdxType = uint32_t, const uint64_t page_size = 4096>
struct PageORAM {
  using Page_ = Page<T, UidType, page_size>;
  using FrontEndType = EM::MemoryServer::NonCachedServerFrontendInstance<
      Page_, EM::Backend::MemServerBackend,
      EM::MemoryServer::EncryptType::ENCRYPT_AND_AUTH_FRESH>;
  using BackendType = typename FrontEndType::BackendType;
  static constexpr uint64_t item_per_page = Page_::item_per_page;

 private:
  std::vector<PageIdxType> posMap;
  struct LinkedNode {
    T data;
    UidType uid;
    UidType nextNodeIdx;
  };
  std::vector<LinkedNode> stash;  // A fairly large stash
  // Each page can have a linked list of nodes
  // Since the size of each node is the same, we avoid using malloc. Instead, we
  // manage the free list by ourself. After eviction, we add the head of the
  // linked list to the free list.
  UidType freeListHead, freeListTail;
  static constexpr UidType listEnd = DUMMY<UidType>();
  PageIdxType numPages;
  std::vector<UidType> pageLists;
  FrontEndType frontend;
  T defaultVal = T();

 public:
  PageORAM(BackendType& _backend = *::EM::Backend::g_DefaultBackend)
      : frontend(_backend) {}

  PageORAM(UidType size,
           BackendType& _backend = *::EM::Backend::g_DefaultBackend)
      : frontend(_backend) {
    SetSize(size);
  }

  void SetSize(UidType size, uint64_t cacheBytes = 0) {
    numPages = divRoundUp(size, item_per_page) * 1.3;
    // std::cout << "size: " << size << " numPages: " << numPages << std::endl;
    frontend.SetSize(numPages);
    posMap.resize(size);
    for (UidType i = 0; i < size; i++) {
      posMap[i] = UniformRandom(numPages - 1);
    }
    stash.reserve(numPages);
    stash.resize(1);
    stash[0].nextNodeIdx = listEnd;
    pageLists.resize(numPages, listEnd);
    freeListHead = 0;
    freeListTail = 0;
  }

  template <class Func>
    requires UpdateFunction<Func, T>
  void Access(UidType address, const Func& accessor) {
    PageIdxType pageIdx = posMap[address];
    Page_ page;
    frontend.Read(pageIdx, page);
    access(address, accessor, pageIdx, page);
    frontend.Write(pageIdx, page);
  }

  /**
   * @brief Read the data at a given address
   *
   * @param address The address to read
   * @param out The output data
   */
  void Read(UidType address, T& out) {
    Access(address, [&](const T& data) { out = data; });
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

  template <typename Reader>
    requires Readable<Reader, T>
  void InitFromReader(Reader& reader,
                      uint64_t cacheBytes = DEFAULT_HEAP_SIZE / 5) {
    // first load everything to the external array
    UidType uid = 0;
    UidType initSize = reader.size();
    for (PageIdxType i = 0; i < numPages; i++) {
      if (reader.eof()) {
        break;
      }
      Page_ page;
      frontend.Read(i, page);
      for (UidType j = 0; j < item_per_page; j++) {
        if (reader.eof()) {
          break;
        }
        page.data[j] = reader.read();
        page.uids[j] = uid++;
      }
      frontend.Write(i, page);
    }
    EM::VirtualVector::VirtualWriter<UidBlock<T, UidType>> overflowWriter(
        initSize, [&](const uint64_t idx, const UidBlock<T, UidType>& block) {
          UidType slot = getFreeSlot();
          stash[slot].data = block.data;
          stash[slot].uid = block.uid;
          addSlotToPageList(posMap[block.uid], slot);
        });
    partitionHelper(0, numPages, cacheBytes, overflowWriter);
  }

  void InitDefault(const T& defaultVal) { this->defaultVal = defaultVal; }

  void PrintMemoryUsage() {
    const double MBinv = 1.0 / 1024 / 1024;
    printf(
        "PageORAM: %lu elements, %lu pages, maximum of %lu elements in stash\n",
        posMap.size(), (uint64_t)numPages, stash.size());
    printf("Raw data size: %f MB\n", posMap.size() * sizeof(T) * MBinv);
    uint64_t stashBytes = stash.size() * sizeof(LinkedNode);
    uint64_t posMapBytes = posMap.size() * sizeof(PageIdxType);
    uint64_t pageListPtrsBytes = pageLists.size() * sizeof(UidType);
    uint64_t freshnessCheckBytes = numPages * sizeof(uint32_t);

    printf("Heap usage: %f MB\n", (stashBytes + posMapBytes +
                                   pageListPtrsBytes + freshnessCheckBytes) *
                                      MBinv);
    printf(
        "(Stash: %f MB, posMap: %f MB, page list heads: %f MB, "
        "freshnessCheck: %f MB)\n",
        stashBytes * MBinv, posMapBytes * MBinv, pageListPtrsBytes * MBinv,
        freshnessCheckBytes * MBinv);
    printf("External memory usage: %f MB\n",
           numPages * sizeof(typename Page_::Encrypted_t) * MBinv);
  }

 private:
  template <class Func>
    requires UpdateFunction<Func, T>
  void access(UidType address, const Func& accessor, PageIdxType pageIdx,
              Page_& page) {
    bool foundInStash = false;
    bool foundInPage = false;
    PageIdxType newPageIdx = UniformRandom(numPages - 1);

    T* dataPtr = nullptr;
    UidType* prevUidPtr = &pageLists[pageIdx];
    UidType currIdx = *prevUidPtr;
    while (currIdx != listEnd) {
      if (stash[currIdx].uid == address) {
        dataPtr = &stash[currIdx].data;
        // remove the node from the page linked list
        *prevUidPtr = stash[currIdx].nextNodeIdx;
        // add the node to the head of the new page linked list
        addSlotToPageList(newPageIdx, currIdx);
        foundInStash = true;
        // break;
      }
      prevUidPtr = &stash[currIdx].nextNodeIdx;
      currIdx = *prevUidPtr;
    }

    for (UidType i = 0; i < item_per_page; i++) {
      if (page.uids[i] == address) {
        dataPtr = &page.data[i];
        page.uids[i] = DUMMY<UidType>();
        foundInPage = true;
        // break;
      }
    }
    if (!foundInStash) {
      UidType freeSlot = getFreeSlot();
      stash[freeSlot].uid = address;
      const T* src = foundInPage ? dataPtr : &defaultVal;

      memcpy(&stash[freeSlot].data, src, sizeof(T));

      dataPtr = &stash[freeSlot].data;
      addSlotToPageList(newPageIdx, freeSlot);
    }

    accessor(*dataPtr);

    // Write back and maintainence

    posMap[address] = newPageIdx;

    uint32_t freeSlotIndices[item_per_page];
    uint32_t numFreeSlots = 0;
    for (UidType i = 0; i < item_per_page; i++) {
      if (page.uids[i] == DUMMY<UidType>()) {
        freeSlotIndices[numFreeSlots++] = i;
      }
    }
    // evict the data of the current page
    UidType prev = listEnd;
    UidType curr = pageLists[pageIdx];
    for (UidType i = 0; i < numFreeSlots; i++) {
      if (curr == listEnd) {
        break;
      }

      uint32_t freeSlotIdx = freeSlotIndices[i];
      page.uids[freeSlotIdx] = stash[curr].uid;
      memcpy(&page.data[freeSlotIdx], &stash[curr].data, sizeof(T));

      prev = curr;
      curr = stash[curr].nextNodeIdx;
    }
    if (prev != listEnd) {
      // add the evicted slots to the free list
      stash[prev].nextNodeIdx = freeListHead;
      freeListHead = pageLists[pageIdx];
    }
    pageLists[pageIdx] = curr;
  }

  template <class OverflowWriter>
  void partitionHelper(PageIdxType beginPageIdx, PageIdxType endPageIdx,
                       uint64_t cacheBytes, OverflowWriter& overflowWriter) {
    PageIdxType pageCounts = endPageIdx - beginPageIdx;
    if (pageCounts <= 1) {
      return;
    }
    PageIdxType maxWay = cacheBytes / sizeof(Page_);
    int levelNeeded = 1;
    for (PageIdxType totalWay = maxWay; totalWay < pageCounts;
         totalWay *= maxWay) {
      levelNeeded++;
    }
    PageIdxType way = ceil(pow(pageCounts, 1.0 / levelNeeded));
    PageIdxType stepSize = cacheBytes / (way * sizeof(Page_));
    if (stepSize == 0) {
      stepSize = 1;
    }
    PageIdxType partitionSize = divRoundUp(pageCounts, way);
    stepSize = std::min(stepSize, partitionSize);
    PageIdxType numSteps = divRoundUp(partitionSize, stepSize);
    {
      std::vector<std::vector<UidBlock<T, UidType>>> partitions(way);
      for (PageIdxType i = 0; i < way; i++) {
        partitions[i].reserve(stepSize * item_per_page * 2);
      }
      std::vector<PageIdxType> stepBeginIndices(way);
      std::vector<PageIdxType> stepEndIndices(way);
      for (PageIdxType i = 0; i < numSteps; i++) {
        for (PageIdxType wayIdx = 0; wayIdx < way; wayIdx++) {
          stepBeginIndices[wayIdx] = std::min(
              beginPageIdx + i * stepSize + wayIdx * partitionSize, endPageIdx);
          stepEndIndices[wayIdx] =
              std::min(std::min(stepBeginIndices[wayIdx] + stepSize,
                                beginPageIdx + (wayIdx + 1) * partitionSize),
                       endPageIdx);
        }
        for (PageIdxType wayIdx = 0; wayIdx < way; wayIdx++) {
          PageIdxType stepBeginPageIdx = stepBeginIndices[wayIdx];
          PageIdxType stepEndPageIdx = stepEndIndices[wayIdx];
          for (PageIdxType j = stepBeginPageIdx; j < stepEndPageIdx; j++) {
            Page_ page;
            frontend.Read(j, page);
            for (uint32_t k = 0; k < item_per_page; k++) {
              UidType uid = page.uids[k];
              if (uid != DUMMY<UidType>()) {
                PageIdxType pageIdx = posMap[uid];
                PageIdxType partitionIdx =
                    (pageIdx - beginPageIdx) / partitionSize;
                UidType partitionMaxSize = (stepEndIndices[partitionIdx] -
                                            stepBeginIndices[partitionIdx]) *
                                           item_per_page;
                // printf("partitionIdx: %u, partitionMaxSize: %lu\n",
                //        partitionIdx, partitionMaxSize);
                if (partitions[partitionIdx].size() < partitionMaxSize) {
                  partitions[partitionIdx].emplace_back(page.data[k], uid);
                } else {
                  overflowWriter.write(UidBlock<T, UidType>(page.data[k], uid));
                }
              }
            }
          }
        }
        for (PageIdxType wayIdx = 0; wayIdx < way; wayIdx++) {
          PageIdxType stepBeginPageIdx = stepBeginIndices[wayIdx];
          PageIdxType stepEndPageIdx = stepEndIndices[wayIdx];
          UidType wayOffset = 0;
          for (PageIdxType j = stepBeginPageIdx; j < stepEndPageIdx; j++) {
            Page_ page = Page_();
            for (uint32_t k = 0; k < item_per_page; k++) {
              if (wayOffset < partitions[wayIdx].size()) {
                page.uids[k] = partitions[wayIdx][wayOffset].uid;
                page.data[k] = partitions[wayIdx][wayOffset].data;
                wayOffset++;
              }
            }
            frontend.Write(j, page);
          }
          Assert(wayOffset == partitions[wayIdx].size());
          partitions[wayIdx].clear();
        }
      }
    }
    for (PageIdxType wayIdx = 0; wayIdx < way; wayIdx++) {
      PageIdxType wayBeginPageIdx = beginPageIdx + wayIdx * partitionSize;
      PageIdxType wayEndPageIdx =
          std::min(wayBeginPageIdx + partitionSize, endPageIdx);
      partitionHelper(wayBeginPageIdx, wayEndPageIdx, cacheBytes,
                      overflowWriter);
    }
  }

  UidType getFreeSlot() {
    if (freeListHead == freeListTail) {  // provision a new slot
      stash.emplace_back();
      freeListTail = stash[freeListTail].nextNodeIdx = stash.size() - 1;
      stash[freeListTail].nextNodeIdx = listEnd;
    }
    UidType slot = freeListHead;
    freeListHead = stash[slot].nextNodeIdx;
    return slot;
  }

  void addSlotToPageList(PageIdxType pageIdx, UidType slot) {
    stash[slot].nextNodeIdx = pageLists[pageIdx];
    pageLists[pageIdx] = slot;
  }

  void printFreeList() {
    printf("Free list: ");
    UidType curr = freeListHead;
    while (curr != listEnd) {
      printf("%lu ", curr);
      curr = stash[curr].nextNodeIdx;
    }
    printf("\n");
  }
};

}  // namespace ODSL