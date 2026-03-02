#ifndef LocalCorePotentialPrimRecKS
#define LocalCorePotentialPrimRecKS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [K|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ks The index of integral in primitive integrals buffer.
/// @param idx_hs The index of integral in primitive integrals buffer.
/// @param idx_is The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_ks(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ks,
                                  const size_t idx_hs,
                                  const size_t idx_is,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecKS */
