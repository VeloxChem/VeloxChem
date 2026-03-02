#ifndef T2CHrrABRecPG_hpp
#define T2CHrrABRecPG_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|G]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pg The index of integral in contracted integrals buffer.
/// @param idx_sg The index of integral in contracted integrals buffer.
/// @param idx_sh The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pg(CSimdArray<double>& cbuffer, 
            const size_t idx_pg,
            const size_t idx_sg,
            const size_t idx_sh,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPG_hpp */

