#ifndef NuclearPotentialGeom010PrimRecHF
#define NuclearPotentialGeom010PrimRecHF

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [H|AG(1)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_hf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_ff The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_ff The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_010_hf(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_hf,
                                        const size_t idx_npot_geom_010_0_ff,
                                        const size_t idx_npot_geom_010_1_ff,
                                        const size_t idx_npot_geom_010_0_gd,
                                        const size_t idx_npot_geom_010_1_gd,
                                        const size_t idx_npot_1_gf,
                                        const size_t idx_npot_geom_010_0_gf,
                                        const size_t idx_npot_geom_010_1_gf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom010PrimRecHF */
