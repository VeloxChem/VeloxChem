#ifndef T2CHrrABRecGD_hpp
#define T2CHrrABRecGD_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [G|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_gd The index of integral in contracted integrals buffer.
/// @param idx_gp The index of integral in contracted integrals buffer.
/// @param idx_hp The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_gd(CSimdArray<double>& cbuffer, 
            const size_t idx_gd,
            const size_t idx_gp,
            const size_t idx_hp,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecGD_hpp */

