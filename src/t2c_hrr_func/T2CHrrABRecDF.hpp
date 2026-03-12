#ifndef T2CHrrABRecDF_hpp
#define T2CHrrABRecDF_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [D|X|F]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_df The index of integral in contracted integrals buffer.
/// @param idx_pf The index of integral in contracted integrals buffer.
/// @param idx_pg The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_df(CSimdArray<double>& cbuffer, 
            const size_t idx_df,
            const size_t idx_pf,
            const size_t idx_pg,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecDF_hpp */

