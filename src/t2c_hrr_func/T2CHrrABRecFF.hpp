#ifndef T2CHrrABRecFF_hpp
#define T2CHrrABRecFF_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [F|X|F]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_ff The index of integral in contracted integrals buffer.
/// @param idx_df The index of integral in contracted integrals buffer.
/// @param idx_dg The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_ff(CSimdArray<double>& cbuffer, 
            const size_t idx_ff,
            const size_t idx_df,
            const size_t idx_dg,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecFF_hpp */

