#ifndef ProjectedCorePotentialPrimRecDGForS
#define ProjectedCorePotentialPrimRecDGForS

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [D|U_l|G]_S integrals with S projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dg_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sg_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_pg_s_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_sg_s_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_pg_s_1_0_0 The index of integral in primitive integrals buffer.
/// @param p The special projector value.
/// @param idx_sg_s_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_pg_s_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_dg_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_dg_s_0_0_0,
                                        const size_t idx_sg_s_0_0_0,
                                        const size_t idx_pg_s_0_0_0,
                                        const size_t idx_sg_s_1_0_0,
                                        const size_t idx_pg_s_1_0_0,
                                        const int p,
                                        const size_t idx_sg_s_0_0_1,
                                        const size_t idx_pg_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecDGForS */
