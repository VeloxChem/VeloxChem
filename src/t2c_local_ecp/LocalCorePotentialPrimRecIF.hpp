#ifndef LocalCorePotentialPrimRecIF
#define LocalCorePotentialPrimRecIF

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [I|U_L|F] integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_if The index of integral in primitive integrals buffer.
/// @param idx_gf The index of integral in primitive integrals buffer.
/// @param idx_hd The index of integral in primitive integrals buffer.
/// @param idx_hf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ra The vector of distances R(RA) = R - A.
/// @param idx_zeta The inverted zeta.
auto
comp_prim_local_core_potential_if(CSimdArray<double>& pbuffer, 
                                  const size_t idx_if,
                                  const size_t idx_gf,
                                  const size_t idx_hd,
                                  const size_t idx_hf,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecIF */
