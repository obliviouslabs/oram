#pragma once

#include "circuit_oram.hpp"
#include "linear_oram.hpp"
#include "path_oram.hpp"

namespace ODSL {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct ORAM {
  using LinearORAM_ = LinearORAM::LinearORAM<T, PositionType, UidType>;
  using ORAM_ = CircuitORAM::ORAM<T, 2, 50, PositionType, UidType>;
  // using ORAM_ = PathORAM::ORAM<T, 5, 64, PositionType, UidType>;
  LinearORAM_* linearOram = NULL;
  ORAM_* pathOram = NULL;
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
      pathOram = new ORAM_(size, cacheBytes);
    }
  }

  static size_t GetMemoryUsage(size_t size, size_t cacheBytes = 1UL << 62) {
    if (size <= linear_oram_threshold) {
      return LinearORAM_::GetMemoryUsage(size);
    } else {
      return ORAM_::GetMemoryUsage(size, cacheBytes);
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

  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
    } else {
      return pathOram->Read(pos, uid, out, newPos);
    }
  }

  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    if (isLinear) {
      return linearOram->Write(uid, in);
    } else {
      return pathOram->Write(uid, in, newPos);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return pathOram->Update(pos, uid, newPos, updateFunc);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc, T& out) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
    } else {
      return pathOram->Update(pos, uid, newPos, updateFunc, out);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return pathOram->Update(pos, uid, newPos, updateFunc, out, updatedUid);
    }
  }

  UidType GetNextUid(bool real = true) {
    UidType res = DUMMY<UidType>();
    obliMove(real, res, nextUid);
    nextUid += (UidType)real;
    return res;
  }
};

}  // namespace ODSL