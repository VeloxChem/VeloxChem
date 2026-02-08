#ifndef T2CHrrABRecHP_hpp
#define T2CHrrABRecHP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [H|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_hp The index of integral in contracted integrals buffer.
/// @param idx_hs The index of integral in contracted integrals buffer.
/// @param idx_is The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_hp(CSimdArray<double>& cbuffer, 
            const size_t idx_hp,
            const size_t idx_hs,
            const size_t idx_is,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecHP_hpp */

