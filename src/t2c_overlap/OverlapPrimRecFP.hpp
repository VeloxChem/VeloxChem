#ifndef OverlapPrimRecFP
#define OverlapPrimRecFP

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [F|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_fp The index of integral in primitive integrals buffer.
/// @param idx_ovl_pp The index of integral in primitive integrals buffer.
/// @param idx_ovl_ds The index of integral in primitive integrals buffer.
/// @param idx_ovl_dp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_fp(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_fp,
                          const size_t              idx_ovl_pp,
                          const size_t              idx_ovl_ds,
                          const size_t              idx_ovl_dp,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecFP */
