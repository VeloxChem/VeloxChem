#ifndef ProjectedCorePotentialPrimRecSFForG
#define ProjectedCorePotentialPrimRecSFForG

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [S|U_l|F]_G integrals with G projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sf_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sp_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_g_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_f_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_sp_g_0_1_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_g_0_1_0 The index of integral in primitive integrals buffer.
/// @param idx_sp_d_0_1_1 The index of integral in primitive integrals buffer.
/// @param idx_sd_d_0_1_1 The index of integral in primitive integrals buffer.
/// @param idx_sd_p_1_1_1 The index of integral in primitive integrals buffer.
/// @param idx_sp_s_1_2_1 The index of integral in primitive integrals buffer.
/// @param idx_sd_s_1_2_1 The index of integral in primitive integrals buffer.
/// @param m The special projector value.
/// @param idx_sp_g_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_sd_g_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_sf_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_sf_g_0_0_0,
                                        const size_t idx_sp_g_0_0_0,
                                        const size_t idx_sd_g_0_0_0,
                                        const size_t idx_sd_f_0_0_1,
                                        const size_t idx_sp_g_0_1_0,
                                        const size_t idx_sd_g_0_1_0,
                                        const size_t idx_sp_d_0_1_1,
                                        const size_t idx_sd_d_0_1_1,
                                        const size_t idx_sd_p_1_1_1,
                                        const size_t idx_sp_s_1_2_1,
                                        const size_t idx_sd_s_1_2_1,
                                        const int m,
                                        const size_t idx_sp_g_0_0_1,
                                        const size_t idx_sd_g_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecSFForG */
