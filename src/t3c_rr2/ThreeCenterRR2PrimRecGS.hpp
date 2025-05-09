#ifndef ThreeCenterRR2PrimRecGS
#define ThreeCenterRR2PrimRecGS

#include "SimdArray.hpp"

namespace t3rr2rec { // t3rr2rec namespace

/// @brief Computes primitive [G|GR.R2(r)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gr_gs The index of integral in primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param idx_g_fs The index of integral in primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param idx_g_gs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r_r2_gs(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_gs,
                  const size_t idx_fs,
                  const size_t idx_g_fs,
                  const size_t idx_gs,
                  const size_t idx_g_gs,
                  const CSimdArray<double>& factors,
                  const size_t idx_rgc,
                  const double a_exp,
                  const double c_exp) -> void;
} // t3rr2rec namespace

#endif /* ThreeCenterRR2PrimRecGS */
