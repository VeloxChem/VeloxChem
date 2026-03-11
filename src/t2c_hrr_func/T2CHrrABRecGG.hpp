#ifndef T2CHrrABRecGG_hpp
#define T2CHrrABRecGG_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [G|X|G]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_gg The index of integral in contracted integrals buffer.
/// @param idx_fg The index of integral in contracted integrals buffer.
/// @param idx_fh The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_gg(CSimdArray<double>& cbuffer, 
            const size_t idx_gg,
            const size_t idx_fg,
            const size_t idx_fh,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecGG_hpp */

