#ifndef ThreeCenterOverlapPrimRecSG
#define ThreeCenterOverlapPrimRecSG

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [S|G(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sg The index of integral in primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_sf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgb The vector of distances R(GB) = G - B.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_sg(CSimdArray<double>& pbuffer, 
                     const size_t idx_sg,
                     const size_t idx_sd,
                     const size_t idx_sf,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecSG */
