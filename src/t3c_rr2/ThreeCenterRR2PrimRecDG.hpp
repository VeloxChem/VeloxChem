#ifndef ThreeCenterRR2PrimRecDG
#define ThreeCenterRR2PrimRecDG

#include "SimdArray.hpp"

namespace t3rr2rec { // t3rr2rec namespace

/// @brief Computes primitive [D|GR.R2(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gr_dg The index of integral in primitive integrals buffer.
/// @param idx_pg The index of integral in primitive integrals buffer.
/// @param idx_g_pg The index of integral in primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param idx_g_df The index of integral in primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_g_dg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r_r2_dg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_dg,
                  const size_t idx_pg,
                  const size_t idx_g_pg,
                  const size_t idx_df,
                  const size_t idx_g_df,
                  const size_t idx_dg,
                  const size_t idx_g_dg,
                  const CSimdArray<double>& factors,
                  const size_t idx_rgc,
                  const double a_exp,
                  const double c_exp) -> void;
} // t3rr2rec namespace

#endif /* ThreeCenterRR2PrimRecDG */
