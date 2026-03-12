#ifndef T2CHrrABRecIP_hpp
#define T2CHrrABRecIP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [I|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_ip The index of integral in contracted integrals buffer.
/// @param idx_is The index of integral in contracted integrals buffer.
/// @param idx_ks The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_ip(CSimdArray<double>& cbuffer, 
            const size_t idx_ip,
            const size_t idx_is,
            const size_t idx_ks,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecIP_hpp */

