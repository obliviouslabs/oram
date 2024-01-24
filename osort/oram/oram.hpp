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
  bool isLinear = false;
  ORAM(PositionType size) {
    isLinear = size <= 100;
    if (isLinear) {
      linearOram = new LinearORAM_(size);
    } else {
      pathOram = new PathORAM_(size);
    }
  }

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
    if (isLinear) {
      linearOram->InitFromReader(reader);
    } else {
      pathOram->InitFromReader(reader, writer);
    }
  }

  size_t size() const {
    if (isLinear) {
      return linearOram->size();
    } else {
      return pathOram->size();
    }
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
    } else {
      return pathOram->Read(pos, uid, out);
    }
  }

  PositionType Write(PositionType pos, const UidType& uid, const T& in) {
    if (isLinear) {
      return linearOram->Write(pos, uid, in);
    } else {
      return pathOram->Write(pos, uid, in);
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
};

}  // namespace ORAM