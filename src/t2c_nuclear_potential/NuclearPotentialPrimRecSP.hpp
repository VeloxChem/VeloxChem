#ifndef NuclearPotentialPrimRecSP
#define NuclearPotentialPrimRecSP

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|A|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param idx_rpc The vector of distances R(PC) = P - C.
auto comp_prim_nuclear_potential_sp(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_sp,
                                    const size_t              idx_npot_0_ss,
                                    const size_t              idx_npot_1_ss,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpb,
                                    const size_t              idx_rpc) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecSP */
