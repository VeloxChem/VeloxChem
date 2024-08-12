#ifndef OverlapPrimRecGH
#define OverlapPrimRecGH

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [G|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_gh The index of integral in primitive integrals buffer.
/// @param idx_ovl_dh The index of integral in primitive integrals buffer.
/// @param idx_ovl_fg The index of integral in primitive integrals buffer.
/// @param idx_ovl_fh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_gh(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_gh,
                          const size_t              idx_ovl_dh,
                          const size_t              idx_ovl_fg,
                          const size_t              idx_ovl_fh,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecGH */
