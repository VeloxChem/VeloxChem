#ifndef NuclearPotentialPrimRecHS
#define NuclearPotentialPrimRecHS

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [H|A|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_hs The index of integral in primitive integrals buffer.
/// @param idx_npot_0_fs The index of integral in primitive integrals buffer.
/// @param idx_npot_1_fs The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gs The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_hs(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_hs,
                                    const size_t              idx_npot_0_fs,
                                    const size_t              idx_npot_1_fs,
                                    const size_t              idx_npot_0_gs,
                                    const size_t              idx_npot_1_gs,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecHS */