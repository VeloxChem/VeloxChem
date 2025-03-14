#ifndef OverlapPrimRecHP
#define OverlapPrimRecHP

#include "SimdArray.hpp"

namespace ovlrec { // ovlrec namespace

/// @brief Computes primitive [H|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_hp The index of integral in primitive integrals buffer.
/// @param idx_ovl_fp The index of integral in primitive integrals buffer.
/// @param idx_ovl_gs The index of integral in primitive integrals buffer.
/// @param idx_ovl_gp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_overlap_hp(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_hp,
                     const size_t idx_ovl_fp,
                     const size_t idx_ovl_gs,
                     const size_t idx_ovl_gp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpa,
                     const double a_exp) -> void;
} // ovlrec namespace

#endif /* OverlapPrimRecHP */
