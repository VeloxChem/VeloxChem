#ifndef NuclearPotentialPrimRecFD
#define NuclearPotentialPrimRecFD

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [F|A|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_fd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_dd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_dd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_fd(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_fd,
                                    const size_t              idx_npot_0_pd,
                                    const size_t              idx_npot_1_pd,
                                    const size_t              idx_npot_0_dp,
                                    const size_t              idx_npot_1_dp,
                                    const size_t              idx_npot_0_dd,
                                    const size_t              idx_npot_1_dd,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecFD */
