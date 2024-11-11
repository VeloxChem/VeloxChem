#ifndef NuclearPotentialGeom010PrimRecGI
#define NuclearPotentialGeom010PrimRecGI

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [G|AG(1)|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_gi The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_di The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_di The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_fh The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_fh The index of integral in primitive integrals buffer.
/// @param idx_npot_1_fi The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_fi The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_fi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_gi(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_gi,
                                             const size_t              idx_npot_geom_010_0_di,
                                             const size_t              idx_npot_geom_010_1_di,
                                             const size_t              idx_npot_geom_010_0_fh,
                                             const size_t              idx_npot_geom_010_1_fh,
                                             const size_t              idx_npot_1_fi,
                                             const size_t              idx_npot_geom_010_0_fi,
                                             const size_t              idx_npot_geom_010_1_fi,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecGI */
