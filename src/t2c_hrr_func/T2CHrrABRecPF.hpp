#ifndef T2CHrrABRecPF_hpp
#define T2CHrrABRecPF_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|F]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pf The index of integral in contracted integrals buffer.
/// @param idx_sf The index of integral in contracted integrals buffer.
/// @param idx_sg The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pf(CSimdArray<double>& cbuffer, 
            const size_t idx_pf,
            const size_t idx_sf,
            const size_t idx_sg,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPF_hpp */

