#ifndef NuclearPotentialGeom010PrimRecHP
#define NuclearPotentialGeom010PrimRecHP

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [H|AG(1)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_hp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_fp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_fp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gs The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gs The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_010_hp(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_hp,
                                        const size_t idx_npot_geom_010_0_fp,
                                        const size_t idx_npot_geom_010_1_fp,
                                        const size_t idx_npot_geom_010_0_gs,
                                        const size_t idx_npot_geom_010_1_gs,
                                        const size_t idx_npot_1_gp,
                                        const size_t idx_npot_geom_010_0_gp,
                                        const size_t idx_npot_geom_010_1_gp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom010PrimRecHP */
