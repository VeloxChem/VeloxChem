#ifndef T2CHrrABRecPP_hpp
#define T2CHrrABRecPP_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|P]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pp The index of integral in contracted integrals buffer.
/// @param idx_sp The index of integral in contracted integrals buffer.
/// @param idx_sd The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pp(CSimdArray<double>& cbuffer, 
            const size_t idx_pp,
            const size_t idx_sp,
            const size_t idx_sd,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPP_hpp */

