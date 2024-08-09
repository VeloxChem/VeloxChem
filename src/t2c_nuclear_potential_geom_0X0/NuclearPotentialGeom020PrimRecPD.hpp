#ifndef NuclearPotentialGeom020PrimRecPD
#define NuclearPotentialGeom020PrimRecPD

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [P|AG(2)|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_020_pd(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_pd,
                                        const size_t idx_npot_geom_020_0_sp,
                                        const size_t idx_npot_geom_020_1_sp,
                                        const size_t idx_npot_geom_010_1_sd,
                                        const size_t idx_npot_geom_020_0_sd,
                                        const size_t idx_npot_geom_020_1_sd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom020PrimRecPD */
