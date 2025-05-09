#ifndef ThreeCenterRR2PrimRecFG
#define ThreeCenterRR2PrimRecFG

#include "SimdArray.hpp"

namespace t3rr2rec { // t3rr2rec namespace

/// @brief Computes primitive [F|GR.R2(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gr_fg The index of integral in primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_g_dg The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_g_ff The index of integral in primitive integrals buffer.
/// @param idx_fg The index of integral in primitive integrals buffer.
/// @param idx_g_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r_r2_fg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_fg,
                  const size_t idx_dg,
                  const size_t idx_g_dg,
                  const size_t idx_ff,
                  const size_t idx_g_ff,
                  const size_t idx_fg,
                  const size_t idx_g_fg,
                  const CSimdArray<double>& factors,
                  const size_t idx_rgc,
                  const double a_exp,
                  const double c_exp) -> void;
} // t3rr2rec namespace

#endif /* ThreeCenterRR2PrimRecFG */
