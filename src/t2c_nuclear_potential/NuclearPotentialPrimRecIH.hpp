#ifndef NuclearPotentialPrimRecIH
#define NuclearPotentialPrimRecIH

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [I|A|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_ih The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gh The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gh The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hg The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hg The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hh The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_ih(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_ih,
                                    const size_t              idx_npot_0_gh,
                                    const size_t              idx_npot_1_gh,
                                    const size_t              idx_npot_0_hg,
                                    const size_t              idx_npot_1_hg,
                                    const size_t              idx_npot_0_hh,
                                    const size_t              idx_npot_1_hh,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecIH */
