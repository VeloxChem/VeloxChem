#ifndef T2CHrrABRecPH_hpp
#define T2CHrrABRecPH_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|H]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_ph The index of integral in contracted integrals buffer.
/// @param idx_sh The index of integral in contracted integrals buffer.
/// @param idx_si The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_ph(CSimdArray<double>& cbuffer, 
            const size_t idx_ph,
            const size_t idx_sh,
            const size_t idx_si,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPH_hpp */

