#ifndef NuclearPotentialGeom010PrimRecSF
#define NuclearPotentialGeom010PrimRecSF

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|AG(1)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_0_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_geom_010_1_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_sf(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_sf,
                                             const size_t              idx_npot_geom_010_0_sp,
                                             const size_t              idx_npot_geom_010_1_sp,
                                             const size_t              idx_npot_1_sd,
                                             const size_t              idx_npot_geom_010_0_sd,
                                             const size_t              idx_npot_geom_010_1_sd,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpb,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecSF */
