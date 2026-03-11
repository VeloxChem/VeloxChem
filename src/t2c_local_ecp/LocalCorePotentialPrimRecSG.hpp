#ifndef LocalCorePotentialPrimRecSG
#define LocalCorePotentialPrimRecSG

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [S|U_L|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sg The index of integral in primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_sf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_sg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sg,
                                  const size_t idx_sd,
                                  const size_t idx_sf,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecSG */
