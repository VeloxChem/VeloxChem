#ifndef ThreeCenterOverlapGradientPrimRecHH
#define ThreeCenterOverlapGradientPrimRecHH

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [H|GX(r)|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_hh The index of integral in primitive integrals buffer.
/// @param idx_gh The index of integral in primitive integrals buffer.
/// @param idx_hg The index of integral in primitive integrals buffer.
/// @param idx_hh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_hh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hh,
                              const size_t idx_gh,
                              const size_t idx_hg,
                              const size_t idx_hh,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecHH */
