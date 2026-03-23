#ifndef LocalCorePotentialPrimRecGS
#define LocalCorePotentialPrimRecGS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [G|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_gs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gs,
                                  const size_t idx_ds,
                                  const size_t idx_fs,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecGS */
