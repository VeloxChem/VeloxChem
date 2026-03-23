#ifndef T2CHrrABRecKP_hpp
#define T2CHrrABRecKP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [K|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_kp The index of integral in contracted integrals buffer.
/// @param idx_ks The index of integral in contracted integrals buffer.
/// @param idx_ls The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_kp(CSimdArray<double>& cbuffer, 
            const size_t idx_kp,
            const size_t idx_ks,
            const size_t idx_ls,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecKP_hpp */

