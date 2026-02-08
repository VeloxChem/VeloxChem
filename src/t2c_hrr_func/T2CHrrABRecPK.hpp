#ifndef T2CHrrABRecPK_hpp
#define T2CHrrABRecPK_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|K]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pk The index of integral in contracted integrals buffer.
/// @param idx_sk The index of integral in contracted integrals buffer.
/// @param idx_sl The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pk(CSimdArray<double>& cbuffer, 
            const size_t idx_pk,
            const size_t idx_sk,
            const size_t idx_sl,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPK_hpp */

