#ifndef ThreeCenterR2PrimRecGF
#define ThreeCenterR2PrimRecGF

#include "SimdArray.hpp"

namespace t3r2rec { // t3r2rec namespace

/// @brief Computes primitive [G|GR2(r)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_gf The index of integral in primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param idx_fd The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_gp The index of integral in primitive integrals buffer.
/// @param idx_gd The index of integral in primitive integrals buffer.
/// @param idx_gf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r2_gf(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gf,
                const size_t idx_df,
                const size_t idx_fd,
                const size_t idx_ff,
                const size_t idx_gp,
                const size_t idx_gd,
                const size_t idx_gf,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void;
} // t3r2rec namespace

#endif /* ThreeCenterR2PrimRecGF */
