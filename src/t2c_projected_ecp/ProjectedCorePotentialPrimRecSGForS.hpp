#ifndef ProjectedCorePotentialPrimRecSGForS
#define ProjectedCorePotentialPrimRecSGForS

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [S|U_l|G]_S integrals with S projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sg_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sf_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sd_s_0_1_0 The index of integral in primitive integrals buffer.
/// @param idx_sf_s_0_1_0 The index of integral in primitive integrals buffer.
/// @param m The special projector value.
/// @param idx_sd_s_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_sf_s_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_sg_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_sg_s_0_0_0,
                                        const size_t idx_sd_s_0_0_0,
                                        const size_t idx_sf_s_0_0_0,
                                        const size_t idx_sd_s_0_1_0,
                                        const size_t idx_sf_s_0_1_0,
                                        const int m,
                                        const size_t idx_sd_s_0_0_1,
                                        const size_t idx_sf_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecSGForS */
