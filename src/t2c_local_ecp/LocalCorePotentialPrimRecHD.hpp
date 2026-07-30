#ifndef LocalCorePotentialPrimRecHD
#define LocalCorePotentialPrimRecHD

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [H|U_L|D] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hd The index of integral in primitive integrals buffer.
/// @param idx_fd The index of integral in primitive integrals buffer.
/// @param idx_gp The index of integral in primitive integrals buffer.
/// @param idx_gd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_hd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hd,
                                  const size_t idx_fd,
                                  const size_t idx_gp,
                                  const size_t idx_gd,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecHD */
