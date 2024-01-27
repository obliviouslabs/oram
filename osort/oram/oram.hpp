#pragma once

#include "linear_oram.hpp"
#include "path_oram.hpp"

namespace ORAM {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct ORAM {
  using LinearORAM_ = LinearORAM::LinearORAM<T, PositionType, UidType>;
  using PathORAM_ = PathORAM::PathORAM<T, 5, 63, PositionType, UidType>;
  LinearORAM_* linearOram = NULL;
  PathORAM_* pathOram = NULL;
  UidType nextUid = 0;
  bool isLinear = false;
  static constexpr PositionType linear_oram_threshold = 100;

  ORAM() {}

  ORAM(PositionType size) { SetSize(size); }

  ORAM(PositionType size, size_t cacheBytes) { SetSize(size, cacheBytes); }

  // TODO add copy constructor

  ~ORAM() {
    if (linearOram) {
      delete linearOram;
      linearOram = NULL;
    }
    if (pathOram) {
      delete pathOram;
      pathOram = NULL;
    }
  }

  template <typename Reader, typename Writer>
  void InitFromReader(Reader& reader, Writer& writer) {
    nextUid = reader.size();
    if (isLinear) {
      linearOram->InitFromReader(reader);
    } else {
      pathOram->InitFromReader(reader, writer);
    }
  }

  PositionType size() const {
    if (isLinear) {
      return linearOram->size();
    } else {
      return pathOram->size();
    }
  }

  void SetSize(PositionType size, size_t cacheBytes = 1UL << 62) {
    if (linearOram || pathOram) {
      throw std::runtime_error("SetSize can only be called on empty oram");
    }
    isLinear = size <= linear_oram_threshold;
    if (isLinear) {
      if (LinearORAM_::GetMemoryUsage(size) > cacheBytes) {
        throw std::runtime_error(
            "LinearORAM_::GetMemoryUsage(size) > cacheBytes");
      }
      linearOram = new LinearORAM_(size);
    } else {
      pathOram = new PathORAM_(size, cacheBytes);
    }
  }

  static size_t GetMemoryUsage(size_t size, size_t cacheBytes = 1UL << 62) {
    if (size <= linear_oram_threshold) {
      return LinearORAM_::GetMemoryUsage(size);
    } else {
      return PathORAM_::GetMemoryUsage(size, cacheBytes);
    }
  }

  size_t GetMemoryUsage() const {
    if (isLinear) {
      Assert(linearOram);
      return linearOram->GetMemoryUsage();
    } else {
      Assert(pathOram);
      return pathOram->GetMemoryUsage();
    }
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
    } else {
      return pathOram->Read(pos, uid, out);
    }
  }

  PositionType Write(const UidType& uid, const T& in) {
    if (isLinear) {
      return linearOram->Write(uid, in);
    } else {
      return pathOram->Write(uid, in);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return pathOram->Update(pos, uid, updateFunc);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
    } else {
      return pathOram->Update(pos, uid, updateFunc, out);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return pathOram->Update(pos, uid, updateFunc, out, updatedUid);
    }
  }

  UidType GetNextUid(bool real = true) {
    UidType res = DUMMY<UidType>();
    obliMove(real, res, nextUid);
    nextUid += (UidType)real;
    return res;
  }
};

}  // namespace ORAM