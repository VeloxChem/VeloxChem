#ifndef NuclearPotentialGeom010PrimRecIF
#define NuclearPotentialGeom010PrimRecIF

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [I|AG(1)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_if The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_hd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_hd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_hf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_hf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_if(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_if,
                                             const size_t              idx_npot_geom_010_0_gf,
                                             const size_t              idx_npot_geom_010_1_gf,
                                             const size_t              idx_npot_geom_010_0_hd,
                                             const size_t              idx_npot_geom_010_1_hd,
                                             const size_t              idx_npot_1_hf,
                                             const size_t              idx_npot_geom_010_0_hf,
                                             const size_t              idx_npot_geom_010_1_hf,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecIF */
