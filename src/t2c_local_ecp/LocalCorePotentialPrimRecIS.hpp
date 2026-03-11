#ifndef LocalCorePotentialPrimRecIS
#define LocalCorePotentialPrimRecIS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [I|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_is The index of integral in primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param idx_hs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_is(CSimdArray<double>& pbuffer, 
                                  const size_t idx_is,
                                  const size_t idx_gs,
                                  const size_t idx_hs,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecIS */
