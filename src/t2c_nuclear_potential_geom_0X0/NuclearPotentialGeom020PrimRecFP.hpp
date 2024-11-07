#ifndef NuclearPotentialGeom020PrimRecFP
#define NuclearPotentialGeom020PrimRecFP

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [F|AG(2)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_020_0_fp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_pp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_pp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_ds The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_ds The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_0_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_020_1_dp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_020_fp(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_020_0_fp,
                                             const size_t              idx_npot_geom_020_0_pp,
                                             const size_t              idx_npot_geom_020_1_pp,
                                             const size_t              idx_npot_geom_020_0_ds,
                                             const size_t              idx_npot_geom_020_1_ds,
                                             const size_t              idx_npot_geom_010_1_dp,
                                             const size_t              idx_npot_geom_020_0_dp,
                                             const size_t              idx_npot_geom_020_1_dp,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom020PrimRecFP */
