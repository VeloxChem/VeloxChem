#ifndef T2CHrrABRecFG_hpp
#define T2CHrrABRecFG_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [F|X|G]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fg The index of integral in contracted integrals buffer.
/// @param idx_dg The index of integral in contracted integrals buffer.
/// @param idx_dh The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_fg(CSimdArray<double>& cbuffer, 
            const size_t idx_fg,
            const size_t idx_dg,
            const size_t idx_dh,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecFG_hpp */

