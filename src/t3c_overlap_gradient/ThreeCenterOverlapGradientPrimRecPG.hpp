#ifndef ThreeCenterOverlapGradientPrimRecPG
#define ThreeCenterOverlapGradientPrimRecPG

#include "SimdArray.hpp"

namespace g3ovlrec { // g3ovlrec namespace

/// @brief Computes primitive [P|GX(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_pg The index of integral in primitive integrals buffer.
/// @param idx_sg The index of integral in primitive integrals buffer.
/// @param idx_pf The index of integral in primitive integrals buffer.
/// @param idx_pg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gradient_pg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_pg,
                              const size_t idx_sg,
                              const size_t idx_pf,
                              const size_t idx_pg,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void;
} // g3ovlrec namespace

#endif /* ThreeCenterOverlapGradientPrimRecPG */
