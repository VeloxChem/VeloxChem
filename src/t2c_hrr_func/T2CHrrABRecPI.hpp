#ifndef T2CHrrABRecPI_hpp
#define T2CHrrABRecPI_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [P|X|I]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pi The index of integral in contracted integrals buffer.
/// @param idx_si The index of integral in contracted integrals buffer.
/// @param idx_sk The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_pi(CSimdArray<double>& cbuffer, 
            const size_t idx_pi,
            const size_t idx_si,
            const size_t idx_sk,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecPI_hpp */

