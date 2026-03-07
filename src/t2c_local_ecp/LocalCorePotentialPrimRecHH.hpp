#ifndef LocalCorePotentialPrimRecHH
#define LocalCorePotentialPrimRecHH

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [H|U_L|H] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hh The index of integral in primitive integrals buffer.
/// @param idx_fh The index of integral in primitive integrals buffer.
/// @param idx_gg The index of integral in primitive integrals buffer.
/// @param idx_gh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_hh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hh,
                                  const size_t idx_fh,
                                  const size_t idx_gg,
                                  const size_t idx_gh,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecHH */
