#ifndef ThreeCenterRR2PrimRecSD
#define ThreeCenterRR2PrimRecSD

#include "SimdArray.hpp"

namespace t3rr2rec { // t3rr2rec namespace

/// @brief Computes primitive [S|GR.R2(r)|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gr_sd The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param idx_g_sp The index of integral in primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_g_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r_r2_sd(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_sd,
                  const size_t idx_sp,
                  const size_t idx_g_sp,
                  const size_t idx_sd,
                  const size_t idx_g_sd,
                  const CSimdArray<double>& factors,
                  const size_t idx_rgc,
                  const double a_exp,
                  const double c_exp) -> void;
} // t3rr2rec namespace

#endif /* ThreeCenterRR2PrimRecSD */
