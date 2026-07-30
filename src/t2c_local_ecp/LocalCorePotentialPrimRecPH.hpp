#ifndef LocalCorePotentialPrimRecPH
#define LocalCorePotentialPrimRecPH

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [P|U_L|H] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ph The index of integral in primitive integrals buffer.
/// @param idx_sg The index of integral in primitive integrals buffer.
/// @param idx_sh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_ph(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ph,
                                  const size_t idx_sg,
                                  const size_t idx_sh,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecPH */
