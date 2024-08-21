#ifndef NuclearPotentialPrimRecII
#define NuclearPotentialPrimRecII

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [I|A|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_ii The index of integral in primitive integrals buffer.
/// @param idx_npot_0_gi The index of integral in primitive integrals buffer.
/// @param idx_npot_1_gi The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hh The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hh The index of integral in primitive integrals buffer.
/// @param idx_npot_0_hi The index of integral in primitive integrals buffer.
/// @param idx_npot_1_hi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_ii(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_ii,
                                    const size_t              idx_npot_0_gi,
                                    const size_t              idx_npot_1_gi,
                                    const size_t              idx_npot_0_hh,
                                    const size_t              idx_npot_1_hh,
                                    const size_t              idx_npot_0_hi,
                                    const size_t              idx_npot_1_hi,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpa,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecII */