#ifndef ThreeCenterOverlapGradientPrimRecPH
#define ThreeCenterOverlapGradientPrimRecPH

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [P|GX(r)|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_ph The index of integral in primitive integrals buffer.
/// @param idx_sh The index of integral in primitive integrals buffer.
/// @param idx_pg The index of integral in primitive integrals buffer.
/// @param idx_ph The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_ph(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ph,
                              const size_t idx_sh,
                              const size_t idx_pg,
                              const size_t idx_ph,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecPH */
