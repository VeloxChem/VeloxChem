#ifndef ThreeCenterR2PrimRecDP
#define ThreeCenterR2PrimRecDP

#include "SimdArray.hpp"

namespace t3r2rec { // t3r2rec namespace

/// @brief Computes primitive [D|GR2(r)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_dp The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param idx_ps The index of integral in primitive integrals buffer.
/// @param idx_pp The index of integral in primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param idx_dp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r2_dp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_dp,
                const size_t idx_sp,
                const size_t idx_ps,
                const size_t idx_pp,
                const size_t idx_ds,
                const size_t idx_dp,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void;
} // t3r2rec namespace

#endif /* ThreeCenterR2PrimRecDP */
