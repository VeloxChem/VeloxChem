#ifndef LocalCorePotentialPrimRecHS
#define LocalCorePotentialPrimRecHS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [H|U_L|S] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hs The index of integral in primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_hs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hs,
                                  const size_t idx_fs,
                                  const size_t idx_gs,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecHS */
