#ifndef OverlapPrimRecPH
#define OverlapPrimRecPH

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [P|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_ph The index of integral in primitive integrals buffer.
/// @param idx_ovl_sg The index of integral in primitive integrals buffer.
/// @param idx_ovl_sh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_ph(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_ph,
                          const size_t              idx_ovl_sg,
                          const size_t              idx_ovl_sh,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecPH */
