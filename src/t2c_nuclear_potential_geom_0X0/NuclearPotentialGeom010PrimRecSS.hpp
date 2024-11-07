#ifndef NuclearPotentialGeom010PrimRecSS_hpp
#define NuclearPotentialGeom010PrimRecSS_hpp

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|AG(1)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_geom_010_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_geom_010_ss(CSimdArray<double>&       pbuffer,
                                             const size_t              idx_npot_geom_010_0_ss,
                                             const size_t              idx_npot_1_ss,
                                             const CSimdArray<double>& factors,
                                             const size_t              idx_rpc,
                                             const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGeom010PrimRecSS_hpp */
