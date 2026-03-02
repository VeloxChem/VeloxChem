#ifndef T2CHrrABRecPD_hpp
#define T2CHrrABRecPD_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pd The index of integral in contracted integrals buffer.
/// @param idx_sd The index of integral in contracted integrals buffer.
/// @param idx_sf The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pd(CSimdArray<double>& cbuffer, 
            const size_t idx_pd,
            const size_t idx_sd,
            const size_t idx_sf,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPD_hpp */

