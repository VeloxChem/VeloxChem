#ifndef NuclearPotentialGeom020PrimRecPS
#define NuclearPotentialGeom020PrimRecPS

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [P|AG(2)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
auto comp_prim_nuclear_potential_geom_020_ps(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_020_0_ps,
                                             const size_t              idx_npot_geom_010_1_ss,
                                             const size_t              idx_npot_geom_020_0_ss,
                                             const size_t              idx_npot_geom_020_1_ss,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom020PrimRecPS */
