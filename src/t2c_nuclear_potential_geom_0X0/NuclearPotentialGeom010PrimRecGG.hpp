#ifndef NuclearPotentialGeom010PrimRecGG
#define NuclearPotentialGeom010PrimRecGG

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [G|AG(1)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_gg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_dg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_dg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_ff The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_ff The index of integral in primitive integrals buffer.
/// @param idx_npot_1_fg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_fg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_gg(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_gg,
                                             const size_t              idx_npot_geom_010_0_dg,
                                             const size_t              idx_npot_geom_010_1_dg,
                                             const size_t              idx_npot_geom_010_0_ff,
                                             const size_t              idx_npot_geom_010_1_ff,
                                             const size_t              idx_npot_1_fg,
                                             const size_t              idx_npot_geom_010_0_fg,
                                             const size_t              idx_npot_geom_010_1_fg,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecGG */
