#ifndef T2CHrrABRecFP_hpp
#define T2CHrrABRecFP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [F|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fp The index of integral in contracted integrals buffer.
/// @param idx_fs The index of integral in contracted integrals buffer.
/// @param idx_gs The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_fp(CSimdArray<double>& cbuffer, 
            const size_t idx_fp,
            const size_t idx_fs,
            const size_t idx_gs,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecFP_hpp */

