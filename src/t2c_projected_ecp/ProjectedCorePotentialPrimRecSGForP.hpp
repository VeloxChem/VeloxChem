#ifndef ProjectedCorePotentialPrimRecSGForP
#define ProjectedCorePotentialPrimRecSGForP

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [S|U_l|G]_P integrals with P projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sg_p_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_p_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sf_p_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sf_s_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_sd_p_0_1_0 The index of integral in primitive integrals buffer.
/// @param idx_sf_p_0_1_0 The index of integral in primitive integrals buffer.
/// @param m The special projector value.
/// @param idx_sd_p_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_sf_p_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_sg_p(CSimdArray<double>& pbuffer, 
                                        const size_t idx_sg_p_0_0_0,
                                        const size_t idx_sd_p_0_0_0,
                                        const size_t idx_sf_p_0_0_0,
                                        const size_t idx_sf_s_0_0_1,
                                        const size_t idx_sd_p_0_1_0,
                                        const size_t idx_sf_p_0_1_0,
                                        const int m,
                                        const size_t idx_sd_p_0_0_1,
                                        const size_t idx_sf_p_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecSGForP */
