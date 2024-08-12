#ifndef NuclearPotentialPrimRecHG
#define NuclearPotentialPrimRecHG

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [H|A|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_hg The index of integral in primitive integrals buffer.
/// @param idx_npot_0_fg The index of integral in primitive integrals buffer.
/// @param idx_npot_1_fg The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gf The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gg The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_hg(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_hg,
                               const size_t idx_npot_0_fg,
                               const size_t idx_npot_1_fg,
                               const size_t idx_npot_0_gf,
                               const size_t idx_npot_1_gf,
                               const size_t idx_npot_0_gg,
                               const size_t idx_npot_1_gg,
                               const CSimdArray<double>& factors,
                               const size_t idx_rpa,
                               const size_t idx_rpc,
                               const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialPrimRecHG */
