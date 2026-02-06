#ifndef LocalCorePotentialPrimRecLS
#define LocalCorePotentialPrimRecLS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [L|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ls The index of integral in primitive integrals buffer.
/// @param idx_is The index of integral in primitive integrals buffer.
/// @param idx_ks The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_ls(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ls,
                                  const size_t idx_is,
                                  const size_t idx_ks,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecLS */
