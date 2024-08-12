#ifndef NuclearPotentialPrimRecID
#define NuclearPotentialPrimRecID

#include "SimdArray.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [I|A|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_id The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_nuclear_potential_id(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_id,
                               const size_t idx_npot_0_gd,
                               const size_t idx_npot_1_gd,
                               const size_t idx_npot_0_hp,
                               const size_t idx_npot_1_hp,
                               const size_t idx_npot_0_hd,
                               const size_t idx_npot_1_hd,
                               const CSimdArray<double>& factors,
                               const size_t idx_rpa,
                               const size_t idx_rpc,
                               const double a_exp) -> void;
} // npotrec namespace

#endif /* NuclearPotentialPrimRecID */
