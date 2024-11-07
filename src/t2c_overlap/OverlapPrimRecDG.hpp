#ifndef OverlapPrimRecDG
#define OverlapPrimRecDG

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [D|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_dg The index of integral in primitive integrals buffer.
/// @param idx_ovl_sg The index of integral in primitive integrals buffer.
/// @param idx_ovl_pf The index of integral in primitive integrals buffer.
/// @param idx_ovl_pg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_dg(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_dg,
                          const size_t              idx_ovl_sg,
                          const size_t              idx_ovl_pf,
                          const size_t              idx_ovl_pg,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecDG */
