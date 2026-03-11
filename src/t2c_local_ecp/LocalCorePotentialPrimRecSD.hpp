#ifndef LocalCorePotentialPrimRecSD
#define LocalCorePotentialPrimRecSD

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [S|U_L|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sd The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_sd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sd,
                                  const size_t idx_ss,
                                  const size_t idx_sp,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecSD */
