#ifndef LocalCorePotentialPrimRecHP
#define LocalCorePotentialPrimRecHP

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [H|U_L|P] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hp The index of integral in primitive integrals buffer.
/// @param idx_fp The index of integral in primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param idx_gp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_hp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hp,
                                  const size_t idx_fp,
                                  const size_t idx_gs,
                                  const size_t idx_gp,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecHP */
