#ifndef ThreeCenterOverlapPrimRecSD
#define ThreeCenterOverlapPrimRecSD

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [S|G(r)|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgb The vector of distances R(GB) = G - B.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_sd(CSimdArray<double>& pbuffer, 
                     const size_t idx_sd,
                     const size_t idx_ss,
                     const size_t idx_sp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecSD */
