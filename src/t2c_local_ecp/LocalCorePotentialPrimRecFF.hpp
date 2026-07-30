#ifndef LocalCorePotentialPrimRecFF
#define LocalCorePotentialPrimRecFF

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [F|U_L|F] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_pf The index of integral in primitive integrals buffer.
/// @param idx_dd The index of integral in primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_ff(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ff,
                                  const size_t idx_pf,
                                  const size_t idx_dd,
                                  const size_t idx_df,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecFF */
