#ifndef LocalCorePotentialPrimRecDI
#define LocalCorePotentialPrimRecDI

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [D|U_L|I] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_di The index of integral in primitive integrals buffer.
/// @param idx_si The index of integral in primitive integrals buffer.
/// @param idx_ph The index of integral in primitive integrals buffer.
/// @param idx_pi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_di(CSimdArray<double>& pbuffer, 
                                  const size_t idx_di,
                                  const size_t idx_si,
                                  const size_t idx_ph,
                                  const size_t idx_pi,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecDI */
