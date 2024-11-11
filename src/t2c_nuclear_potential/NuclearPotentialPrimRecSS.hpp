#ifndef NuclearPotentialPrimRecSS
#define NuclearPotentialPrimRecSS

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|A|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_ss The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss  The index of integral in primitive integrals buffer.
/// @param bf_data The Boys function data.
/// @param idx_vals The primary row index of values in Boys function data.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_ss(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_ss,
                                    const size_t              idx_ovl_ss,
                                    const CSimdArray<double>& bf_data,
                                    const size_t              idx_vals,
                                    CSimdArray<double>&       factors,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecSS */
