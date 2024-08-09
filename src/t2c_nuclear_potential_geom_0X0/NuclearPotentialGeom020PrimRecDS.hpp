#ifndef NuclearPotentialGeom020PrimRecDS
#define NuclearPotentialGeom020PrimRecDS

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [D|AG(2)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_ds The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_ps The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_geom_020_ds(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_ds,
                                        const size_t idx_npot_geom_020_0_ss,
                                        const size_t idx_npot_geom_020_1_ss,
                                        const size_t idx_npot_geom_010_1_ps,
                                        const size_t idx_npot_geom_020_0_ps,
                                        const size_t idx_npot_geom_020_1_ps,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGeom020PrimRecDS */
