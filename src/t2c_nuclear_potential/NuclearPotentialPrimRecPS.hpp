#ifndef NuclearPotentialPrimRecPS
#define NuclearPotentialPrimRecPS

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [P|A|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
auto comp_prim_nuclear_potential_ps(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_ps,
                                    const size_t              idx_npot_0_ss,
                                    const size_t              idx_npot_1_ss,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecPS */
