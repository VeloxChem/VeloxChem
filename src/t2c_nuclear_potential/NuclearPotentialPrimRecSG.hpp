#ifndef NuclearPotentialPrimRecSG
#define NuclearPotentialPrimRecSG

#include "SimdArray.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|A|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_npot_0_sg The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_nuclear_potential_sg(CSimdArray<double>&       pbuffer,
                                    const size_t              idx_npot_0_sg,
                                    const size_t              idx_npot_0_sd,
                                    const size_t              idx_npot_1_sd,
                                    const size_t              idx_npot_0_sf,
                                    const size_t              idx_npot_1_sf,
                                    const CSimdArray<double>& factors,
                                    const size_t              idx_rpb,
                                    const size_t              idx_rpc,
                                    const double              a_exp) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialPrimRecSG */
