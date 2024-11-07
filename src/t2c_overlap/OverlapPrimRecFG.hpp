#ifndef OverlapPrimRecFG
#define OverlapPrimRecFG

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [F|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_fg The index of integral in primitive integrals buffer.
/// @param idx_ovl_pg The index of integral in primitive integrals buffer.
/// @param idx_ovl_df The index of integral in primitive integrals buffer.
/// @param idx_ovl_dg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_fg(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_fg,
                          const size_t              idx_ovl_pg,
                          const size_t              idx_ovl_df,
                          const size_t              idx_ovl_dg,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecFG */
