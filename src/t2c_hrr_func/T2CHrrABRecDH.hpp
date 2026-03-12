#ifndef T2CHrrABRecDH_hpp
#define T2CHrrABRecDH_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|H]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dh The index of integral in contracted integrals buffer.
/// @param idx_ph The index of integral in contracted integrals buffer.
/// @param idx_pi The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_dh(CSimdArray<double>& cbuffer, 
            const size_t idx_dh,
            const size_t idx_ph,
            const size_t idx_pi,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDH_hpp */

