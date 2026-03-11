#ifndef T2CHrrABRecFH_hpp
#define T2CHrrABRecFH_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [F|X|H]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fh The index of integral in contracted integrals buffer.
/// @param idx_dh The index of integral in contracted integrals buffer.
/// @param idx_di The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_fh(CSimdArray<double>& cbuffer, 
            const size_t idx_fh,
            const size_t idx_dh,
            const size_t idx_di,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecFH_hpp */

