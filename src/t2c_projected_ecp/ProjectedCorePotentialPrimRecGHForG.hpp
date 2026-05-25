#ifndef ProjectedCorePotentialPrimRecGHForG
#define ProjectedCorePotentialPrimRecGHForG

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [G|U_l|H]_G integrals with G projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gh_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_dh_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_fh_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_fg_f_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_fh_f_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_dh_g_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_fh_g_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_dh_d_1_0_1 The index of integral in primitive integrals buffer.
/// @param idx_fh_d_1_0_1 The index of integral in primitive integrals buffer.
/// @param idx_fg_p_1_1_1 The index of integral in primitive integrals buffer.
/// @param idx_fh_p_1_1_1 The index of integral in primitive integrals buffer.
/// @param idx_dh_s_2_1_1 The index of integral in primitive integrals buffer.
/// @param idx_fh_s_2_1_1 The index of integral in primitive integrals buffer.
/// @param p The special projector value.
/// @param idx_dh_g_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_fh_g_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_gh_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gh_g_0_0_0,
                                        const size_t idx_dh_g_0_0_0,
                                        const size_t idx_fh_g_0_0_0,
                                        const size_t idx_fg_f_0_0_1,
                                        const size_t idx_fh_f_0_0_1,
                                        const size_t idx_dh_g_1_0_0,
                                        const size_t idx_fh_g_1_0_0,
                                        const size_t idx_dh_d_1_0_1,
                                        const size_t idx_fh_d_1_0_1,
                                        const size_t idx_fg_p_1_1_1,
                                        const size_t idx_fh_p_1_1_1,
                                        const size_t idx_dh_s_2_1_1,
                                        const size_t idx_fh_s_2_1_1,
                                        const int p,
                                        const size_t idx_dh_g_0_0_1,
                                        const size_t idx_fh_g_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecGHForG */
