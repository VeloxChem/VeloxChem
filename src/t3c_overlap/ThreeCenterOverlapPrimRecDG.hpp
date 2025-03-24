#ifndef ThreeCenterOverlapPrimRecDG
#define ThreeCenterOverlapPrimRecDG

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [D|G(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_sg The index of integral in primitive integrals buffer.
/// @param idx_pf The index of integral in primitive integrals buffer.
/// @param idx_pg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_dg(CSimdArray<double>& pbuffer, 
                     const size_t idx_dg,
                     const size_t idx_sg,
                     const size_t idx_pf,
                     const size_t idx_pg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecDG */
