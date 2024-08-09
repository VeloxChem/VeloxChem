#ifndef OverlapPrimRecSP
#define OverlapPrimRecSP

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [S|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_sp The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
auto comp_prim_overlap_sp(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_sp,
                          const size_t              idx_ovl_ss,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpb) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecSP */
