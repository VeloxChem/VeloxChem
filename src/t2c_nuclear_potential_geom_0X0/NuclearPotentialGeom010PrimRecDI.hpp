#ifndef NuclearPotentialGeom010PrimRecDI
#define NuclearPotentialGeom010PrimRecDI

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [D|AG(1)|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_di The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_si The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_si The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_ph The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_ph The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pi The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_pi The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_pi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_di(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_di,
                                             const size_t              idx_npot_geom_010_0_si,
                                             const size_t              idx_npot_geom_010_1_si,
                                             const size_t              idx_npot_geom_010_0_ph,
                                             const size_t              idx_npot_geom_010_1_ph,
                                             const size_t              idx_npot_1_pi,
                                             const size_t              idx_npot_geom_010_0_pi,
                                             const size_t              idx_npot_geom_010_1_pi,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecDI */
