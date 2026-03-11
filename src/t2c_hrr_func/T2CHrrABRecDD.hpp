#ifndef T2CHrrABRecDD_hpp
#define T2CHrrABRecDD_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dd The index of integral in contracted integrals buffer.
/// @param idx_pd The index of integral in contracted integrals buffer.
/// @param idx_pf The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_dd(CSimdArray<double>& cbuffer, 
            const size_t idx_dd,
            const size_t idx_pd,
            const size_t idx_pf,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDD_hpp */

