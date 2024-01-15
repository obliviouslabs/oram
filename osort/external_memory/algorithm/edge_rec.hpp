#pragma once
/// @brief EdgeRec is a bitset that records the edges in a graph
/// @tparam Bits the type of the bitset, default is uint64_t
template <typename Bits = uint64_t>
struct EdgeRec {
  Bits edges;
  uint16_t k;

 public:
  EdgeRec(uint16_t k) : k(k) {  // k is number of vertices
    static_assert(std::is_same<Bits, uint64_t>::value,
                  "unrecognized type for recording edge");
    edges = 0;
  }
  // map vertex to edge offset in Edge Rec
  // it's ordered by 0-0, 1-0, 1-1, 2-0, 2-1, ...
  INLINE static constexpr uint16_t getEdgeOffset(uint16_t v0, uint16_t v1) {
    return (v0 << 3) | v1;
  }

  INLINE static constexpr Bits getEdgesMask(uint16_t v0, uint16_t v1) {
    if constexpr (std::is_same<Bits, uint64_t>::value) {
      return (1UL << getEdgeOffset(v0, v1)) | (1UL << getEdgeOffset(v1, v0));
    }
    return 0;
  }

  INLINE void flipEdge(Bits mask) {
    if constexpr (std::is_same<Bits, uint64_t>::value) {
      edges ^= mask;
    }
  }

  INLINE void flipEdge(uint16_t v0, uint16_t v1) {
    flipEdge(getEdgesMask(v0, v1));
  }

  INLINE bool retrieveEdge(uint16_t edgeOffset) const {
    if constexpr (std::is_same<Bits, uint64_t>::value) {
      return (edges >> edgeOffset) & 1UL;
    }
  }

  INLINE bool retrieveEdge(uint16_t v0, uint16_t v1) const {
    uint16_t offset = getEdgeOffset(v0, v1);
    return retrieveEdge(offset);
  }

  // retrieve the direction of v0->v1 and flip the direction
  INLINE bool retrieveAndFlipEdge(uint16_t v0, uint16_t v1) {
    Bits mask = getEdgesMask(v0, v1);
    Bits v0_to_v1_mask = 1UL << getEdgeOffset(v0, v1);
    Bits v1_to_v0_mask = 1UL << getEdgeOffset(v1, v0);
    edges ^= v0_to_v1_mask | v1_to_v0_mask;
    return !(edges & v0_to_v1_mask);
  }

  void print() {
    for (uint16_t i = 0; i < k; ++i) {
      for (uint16_t j = i + 1; j < k; ++j) {
        if (retrieveEdge(i, j)) {
          printf("%d - %d\n", i, j);
        }
      }
    }
  }

  void printPath(const EdgeRec& path) {
    Assert(k == path.k);
    for (uint16_t i = 0; i < k; ++i) {
      for (uint16_t j = i + 1; j < k; ++j) {
        if (retrieveEdge(i, j)) {
          if (path.retrieveEdge(i, j)) {
            printf("%d -> %d\n", i, j);
          } else {
            printf("%d -> %d\n", j, i);
          }
        }
      }
    }
  }

  INLINE EdgeRec EulerPath(uint64_t numEdge) {
    numEdge = std::min(numEdge, (uint64_t)(k * (k - 1) / 2));
    // an edge is set to 1 if it's direction is from small vertex to larger
    EdgeRec path(k);
    // for edges appearing even number of times, set direction according to
    // their order
    path.edges = (~edges) & 0xFF7F3F1F0F070301UL;
    static constexpr Bits simpleMask = ~0x8040201008040201UL;
    // filtering out edges pointing to itself
    edges &= simpleMask;
    uint16_t v8 = 0;
    for (uint16_t r = 0; r != numEdge; ++r) {
      // clear bits for less than v0
      uint64_t choices = edges & ((-1UL) << v8);
      uint16_t trailingZeros = std::countr_zero(choices);
      // case 1: there's an edge v0->next
      // case 2: next is a vertex connected to the smallest uncovered vertex
      // case 3: all is done and next = 0
      uint16_t next = trailingZeros & 7;
      uint16_t v0 = trailingZeros >> 3;
      Bits v0_to_next_mask = 1UL << trailingZeros;
      v8 = next << 3;
      Bits next_to_v0_mask = 1UL << (v8 | v0);
      // set this edge as done (could be that it's already done in case 2)
      edges &= ~(v0_to_next_mask | next_to_v0_mask);
      // set direction v0->next if it's in case 1
      path.edges |= v0_to_next_mask;
    }
    return path;
  }
};