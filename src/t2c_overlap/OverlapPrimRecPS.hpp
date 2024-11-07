#ifndef OverlapPrimRecPS
#define OverlapPrimRecPS

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [P|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_ps The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
auto comp_prim_overlap_ps(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_ps,
                          const size_t              idx_ovl_ss,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecPS */
