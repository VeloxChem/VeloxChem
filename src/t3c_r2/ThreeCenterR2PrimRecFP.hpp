#ifndef ThreeCenterR2PrimRecFP
#define ThreeCenterR2PrimRecFP

#include "SimdArray.hpp"

namespace t3r2rec { // t3r2rec namespace

/// @brief Computes primitive [F|GR2(r)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_g_fp The index of integral in primitive integrals buffer.
/// @param idx_pp The index of integral in primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param idx_dp The index of integral in primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param idx_fp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgc The vector of distances R(GC) = G - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_r2_fp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_fp,
                const size_t idx_pp,
                const size_t idx_ds,
                const size_t idx_dp,
                const size_t idx_fs,
                const size_t idx_fp,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void;
} // t3r2rec namespace

#endif /* ThreeCenterR2PrimRecFP */
