#ifndef ThreeCenterOverlapGradientPrimRecHP
#define ThreeCenterOverlapGradientPrimRecHP

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [H|GX(r)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_hp The index of integral in primitive integrals buffer.
/// @param idx_gp The index of integral in primitive integrals buffer.
/// @param idx_hs The index of integral in primitive integrals buffer.
/// @param idx_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_hp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hp,
                              const size_t idx_gp,
                              const size_t idx_hs,
                              const size_t idx_hp,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecHP */
