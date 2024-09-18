#ifndef NuclearPotentialGeom010PrimRecDG
#define NuclearPotentialGeom010PrimRecDG

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [D|AG(1)|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_dg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_sg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_sg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_pg The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_pg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_dg(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_dg,
                                             const size_t              idx_npot_geom_010_0_sg,
                                             const size_t              idx_npot_geom_010_1_sg,
                                             const size_t              idx_npot_geom_010_0_pf,
                                             const size_t              idx_npot_geom_010_1_pf,
                                             const size_t              idx_npot_1_pg,
                                             const size_t              idx_npot_geom_010_0_pg,
                                             const size_t              idx_npot_geom_010_1_pg,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecDG */