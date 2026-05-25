#ifndef ProjectedCorePotentialPrimRecHGForD
#define ProjectedCorePotentialPrimRecHGForD

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [H|U_l|G]_D integrals with D projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hg_d_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_fg_d_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_gg_d_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_gf_p_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_gg_p_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_fg_d_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_gg_d_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_fg_s_1_0_1 The index of integral in primitive integrals buffer.
/// @param idx_gg_s_1_0_1 The index of integral in primitive integrals buffer.
/// @param p The special projector value.
/// @param idx_fg_d_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_gg_d_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_hg_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_hg_d_0_0_0,
                                        const size_t idx_fg_d_0_0_0,
                                        const size_t idx_gg_d_0_0_0,
                                        const size_t idx_gf_p_0_0_1,
                                        const size_t idx_gg_p_0_0_1,
                                        const size_t idx_fg_d_1_0_0,
                                        const size_t idx_gg_d_1_0_0,
                                        const size_t idx_fg_s_1_0_1,
                                        const size_t idx_gg_s_1_0_1,
                                        const int p,
                                        const size_t idx_fg_d_0_0_1,
                                        const size_t idx_gg_d_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecHGForD */
