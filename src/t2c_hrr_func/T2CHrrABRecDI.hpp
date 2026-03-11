#ifndef T2CHrrABRecDI_hpp
#define T2CHrrABRecDI_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|I]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_di The index of integral in contracted integrals buffer.
/// @param idx_pi The index of integral in contracted integrals buffer.
/// @param idx_pk The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_di(CSimdArray<double>& cbuffer, 
            const size_t idx_di,
            const size_t idx_pi,
            const size_t idx_pk,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDI_hpp */

