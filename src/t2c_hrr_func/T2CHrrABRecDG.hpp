#ifndef T2CHrrABRecDG_hpp
#define T2CHrrABRecDG_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|G]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dg The index of integral in contracted integrals buffer.
/// @param idx_pg The index of integral in contracted integrals buffer.
/// @param idx_ph The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_dg(CSimdArray<double>& cbuffer, 
            const size_t idx_dg,
            const size_t idx_pg,
            const size_t idx_ph,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDG_hpp */

