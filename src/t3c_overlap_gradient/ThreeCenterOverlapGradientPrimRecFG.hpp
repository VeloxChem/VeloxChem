#ifndef ThreeCenterOverlapGradientPrimRecFG
#define ThreeCenterOverlapGradientPrimRecFG

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [F|GX(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_fg The index of integral in primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_fg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fg,
                              const size_t idx_dg,
                              const size_t idx_ff,
                              const size_t idx_fg,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecFG */
