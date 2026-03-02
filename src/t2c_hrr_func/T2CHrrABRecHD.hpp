#ifndef T2CHrrABRecHD_hpp
#define T2CHrrABRecHD_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [H|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_hd The index of integral in contracted integrals buffer.
/// @param idx_hp The index of integral in contracted integrals buffer.
/// @param idx_ip The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_hd(CSimdArray<double>& cbuffer, 
            const size_t idx_hd,
            const size_t idx_hp,
            const size_t idx_ip,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecHD_hpp */

