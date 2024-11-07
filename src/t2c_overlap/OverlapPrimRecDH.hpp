#ifndef OverlapPrimRecDH
#define OverlapPrimRecDH

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [D|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_dh The index of integral in primitive integrals buffer.
/// @param idx_ovl_sh The index of integral in primitive integrals buffer.
/// @param idx_ovl_pg The index of integral in primitive integrals buffer.
/// @param idx_ovl_ph The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_dh(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_dh,
                          const size_t              idx_ovl_sh,
                          const size_t              idx_ovl_pg,
                          const size_t              idx_ovl_ph,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecDH */
