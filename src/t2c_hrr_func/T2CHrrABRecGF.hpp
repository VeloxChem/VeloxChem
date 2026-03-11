#ifndef T2CHrrABRecGF_hpp
#define T2CHrrABRecGF_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [G|X|F]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_gf The index of integral in contracted integrals buffer.
/// @param idx_gd The index of integral in contracted integrals buffer.
/// @param idx_hd The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_gf(CSimdArray<double>& cbuffer, 
            const size_t idx_gf,
            const size_t idx_gd,
            const size_t idx_hd,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecGF_hpp */

