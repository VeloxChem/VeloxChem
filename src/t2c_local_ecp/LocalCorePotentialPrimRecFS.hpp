#ifndef LocalCorePotentialPrimRecFS
#define LocalCorePotentialPrimRecFS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [F|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param idx_ps The index of integral in primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_fs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fs,
                                  const size_t idx_ps,
                                  const size_t idx_ds,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecFS */
