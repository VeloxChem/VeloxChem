#ifndef NuclearPotentialGeom020PrimRecDF
#define NuclearPotentialGeom020PrimRecDF

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [D|AG(2)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_df The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_pf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_020_df(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_df,
                                        const size_t idx_npot_geom_020_0_sf,
                                        const size_t idx_npot_geom_020_1_sf,
                                        const size_t idx_npot_geom_020_0_pd,
                                        const size_t idx_npot_geom_020_1_pd,
                                        const size_t idx_npot_geom_010_1_pf,
                                        const size_t idx_npot_geom_020_0_pf,
                                        const size_t idx_npot_geom_020_1_pf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom020PrimRecDF */
