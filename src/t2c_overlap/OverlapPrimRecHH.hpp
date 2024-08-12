#ifndef OverlapPrimRecHH
#define OverlapPrimRecHH

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [H|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_hh The index of integral in primitive integrals buffer.
/// @param idx_ovl_fh The index of integral in primitive integrals buffer.
/// @param idx_ovl_gg The index of integral in primitive integrals buffer.
/// @param idx_ovl_gh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_hh(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_hh,
                          const size_t              idx_ovl_fh,
                          const size_t              idx_ovl_gg,
                          const size_t              idx_ovl_gh,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecHH */
