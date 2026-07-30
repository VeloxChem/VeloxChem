#ifndef LocalCorePotentialPrimRecDF
#define LocalCorePotentialPrimRecDF

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [D|U_L|F] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param idx_sf The index of integral in primitive integrals buffer.
/// @param idx_pd The index of integral in primitive integrals buffer.
/// @param idx_pf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_df(CSimdArray<double>& pbuffer, 
                                  const size_t idx_df,
                                  const size_t idx_sf,
                                  const size_t idx_pd,
                                  const size_t idx_pf,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecDF */
