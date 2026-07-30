#ifndef LocalCorePotentialPrimRecDS
#define LocalCorePotentialPrimRecDS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [D|U_L|S] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param idx_ps The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_ds(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ds,
                                  const size_t idx_ss,
                                  const size_t idx_ps,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecDS */
