#ifndef NuclearPotentialGeom020PrimRecPG
#define NuclearPotentialGeom020PrimRecPG

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [P|AG(2)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_pg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_sg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_sg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_sg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_020_pg(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_pg,
                                        const size_t idx_npot_geom_020_0_sf,
                                        const size_t idx_npot_geom_020_1_sf,
                                        const size_t idx_npot_geom_010_1_sg,
                                        const size_t idx_npot_geom_020_0_sg,
                                        const size_t idx_npot_geom_020_1_sg,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom020PrimRecPG */
