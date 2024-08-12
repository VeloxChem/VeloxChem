#ifndef OverlapPrimRecID
#define OverlapPrimRecID

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [I|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_id The index of integral in primitive integrals buffer.
/// @param idx_ovl_gd The index of integral in primitive integrals buffer.
/// @param idx_ovl_hp The index of integral in primitive integrals buffer.
/// @param idx_ovl_hd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_id(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_id,
                          const size_t              idx_ovl_gd,
                          const size_t              idx_ovl_hp,
                          const size_t              idx_ovl_hd,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecID */
