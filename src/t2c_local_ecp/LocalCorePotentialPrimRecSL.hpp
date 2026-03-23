#ifndef LocalCorePotentialPrimRecSL
#define LocalCorePotentialPrimRecSL

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [S|U_L|L]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sl The index of integral in primitive integrals buffer.
/// @param idx_si The index of integral in primitive integrals buffer.
/// @param idx_sk The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
auto
comp_prim_local_core_potential_sl(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sl,
                                  const size_t idx_si,
                                  const size_t idx_sk,
                                  const CSimdArray<double>& factors) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecSL */
