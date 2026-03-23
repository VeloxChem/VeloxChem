#ifndef T2CHrrABRecHF_hpp
#define T2CHrrABRecHF_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [H|X|F]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_hf The index of integral in contracted integrals buffer.
/// @param idx_hd The index of integral in contracted integrals buffer.
/// @param idx_id The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_hf(CSimdArray<double>& cbuffer, 
            const size_t idx_hf,
            const size_t idx_hd,
            const size_t idx_id,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecHF_hpp */

