#ifndef NuclearPotentialPrimRecIG
#define NuclearPotentialPrimRecIG

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [I|A|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_ig The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gg The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gg The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hf The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hg The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_ig(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_ig,
                                    const size_t              idx_npot_0_gg,
                                    const size_t              idx_npot_1_gg,
                                    const size_t              idx_npot_0_hf,
                                    const size_t              idx_npot_1_hf,
                                    const size_t              idx_npot_0_hg,
                                    const size_t              idx_npot_1_hg,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecIG */
