#ifndef LocalCorePotentialPrimRecPD
#define LocalCorePotentialPrimRecPD

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [P|U_L|D] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_pd The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_pd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pd,
                                  const size_t idx_sp,
                                  const size_t idx_sd,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecPD */
