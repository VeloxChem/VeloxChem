#ifndef ThreeCenterR2PrimRecGG
#define ThreeCenterR2PrimRecGG

#include "SimdArray.hpp"

namespace t3r2rec { // t3r2rec namespace

/// @brief Computes primitive [G|GR2(r)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_gg The index of integral in primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_fg The index of integral in primitive integrals buffer.
/// @param idx_gd The index of integral in primitive integrals buffer.
/// @param idx_gf The index of integral in primitive integrals buffer.
/// @param idx_gg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r2_gg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gg,
                const size_t idx_dg,
                const size_t idx_ff,
                const size_t idx_fg,
                const size_t idx_gd,
                const size_t idx_gf,
                const size_t idx_gg,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void;
} // t3r2rec namespace

#endif /* ThreeCenterR2PrimRecGG */
