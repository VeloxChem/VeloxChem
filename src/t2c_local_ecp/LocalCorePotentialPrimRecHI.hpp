#ifndef LocalCorePotentialPrimRecHI
#define LocalCorePotentialPrimRecHI

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [H|U_L|I] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hi The index of integral in primitive integrals buffer.
/// @param idx_fi The index of integral in primitive integrals buffer.
/// @param idx_gh The index of integral in primitive integrals buffer.
/// @param idx_gi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_hi(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hi,
                                  const size_t idx_fi,
                                  const size_t idx_gh,
                                  const size_t idx_gi,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecHI */
