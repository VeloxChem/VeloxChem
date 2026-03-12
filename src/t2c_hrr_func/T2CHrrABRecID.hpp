#ifndef T2CHrrABRecID_hpp
#define T2CHrrABRecID_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [I|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_id The index of integral in contracted integrals buffer.
/// @param idx_ip The index of integral in contracted integrals buffer.
/// @param idx_kp The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_id(CSimdArray<double>& cbuffer, 
            const size_t idx_id,
            const size_t idx_ip,
            const size_t idx_kp,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecID_hpp */

