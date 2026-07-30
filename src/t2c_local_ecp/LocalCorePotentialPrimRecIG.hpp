#ifndef LocalCorePotentialPrimRecIG
#define LocalCorePotentialPrimRecIG

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [I|U_L|G] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ig The index of integral in primitive integrals buffer.
/// @param idx_gg The index of integral in primitive integrals buffer.
/// @param idx_hf The index of integral in primitive integrals buffer.
/// @param idx_hg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_ig(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ig,
                                  const size_t idx_gg,
                                  const size_t idx_hf,
                                  const size_t idx_hg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecIG */
