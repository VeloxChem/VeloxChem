#ifndef ThreeCenterOverlapGradientPrimRecHF
#define ThreeCenterOverlapGradientPrimRecHF

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [H|GX(r)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_hf The index of integral in primitive integrals buffer.
/// @param idx_gf The index of integral in primitive integrals buffer.
/// @param idx_hd The index of integral in primitive integrals buffer.
/// @param idx_hf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_hf(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hf,
                              const size_t idx_gf,
                              const size_t idx_hd,
                              const size_t idx_hf,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecHF */
