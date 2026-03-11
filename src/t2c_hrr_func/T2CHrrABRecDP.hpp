#ifndef T2CHrrABRecDP_hpp
#define T2CHrrABRecDP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dp The index of integral in contracted integrals buffer.
/// @param idx_ds The index of integral in contracted integrals buffer.
/// @param idx_fs The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_dp(CSimdArray<double>& cbuffer, 
            const size_t idx_dp,
            const size_t idx_ds,
            const size_t idx_fs,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDP_hpp */

