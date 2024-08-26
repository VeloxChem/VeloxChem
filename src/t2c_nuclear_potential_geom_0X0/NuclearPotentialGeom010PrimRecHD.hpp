#ifndef NuclearPotentialGeom010PrimRecHD
#define NuclearPotentialGeom010PrimRecHD

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [H|AG(1)|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_hd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_fd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_fd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_gd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_hd(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_hd,
                                             const size_t              idx_npot_geom_010_0_fd,
                                             const size_t              idx_npot_geom_010_1_fd,
                                             const size_t              idx_npot_geom_010_0_gp,
                                             const size_t              idx_npot_geom_010_1_gp,
                                             const size_t              idx_npot_1_gd,
                                             const size_t              idx_npot_geom_010_0_gd,
                                             const size_t              idx_npot_geom_010_1_gd,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpa,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecHD */
