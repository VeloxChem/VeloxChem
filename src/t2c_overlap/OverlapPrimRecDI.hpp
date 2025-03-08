#ifndef OverlapPrimRecDI
#define OverlapPrimRecDI

#include "SimdArray.hpp"

namespace ovlrec { // ovlrec namespace

/// @brief Computes primitive [D|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_di The index of integral in primitive integrals buffer.
/// @param idx_ovl_si The index of integral in primitive integrals buffer.
/// @param idx_ovl_ph The index of integral in primitive integrals buffer.
/// @param idx_ovl_pi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_overlap_di(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_di,
                     const size_t idx_ovl_si,
                     const size_t idx_ovl_ph,
                     const size_t idx_ovl_pi,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpa,
                     const double a_exp) -> void;
} // ovlrec namespace

#endif /* OverlapPrimRecDI */
