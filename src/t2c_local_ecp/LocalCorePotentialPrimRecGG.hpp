#ifndef LocalCorePotentialPrimRecGG
#define LocalCorePotentialPrimRecGG

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [G|U_L|G] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gg The index of integral in primitive integrals buffer.
/// @param idx_dg The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_gg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gg,
                                  const size_t idx_dg,
                                  const size_t idx_ff,
                                  const size_t idx_fg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecGG */
