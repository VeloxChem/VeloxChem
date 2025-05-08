#ifndef ThreeCenterR2PrimRecSF
#define ThreeCenterR2PrimRecSF

#include "SimdArray.hpp"

namespace t3r2rec { // t3r2rec namespace

/// @brief Computes primitive [S|GR2(r)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_sf The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_sf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r2_sf(CSimdArray<double>& pbuffer, 
                const size_t idx_g_sf,
                const size_t idx_sp,
                const size_t idx_sd,
                const size_t idx_sf,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void;
} // t3r2rec namespace

#endif /* ThreeCenterR2PrimRecSF */
