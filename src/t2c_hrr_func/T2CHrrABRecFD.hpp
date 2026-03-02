#ifndef T2CHrrABRecFD_hpp
#define T2CHrrABRecFD_hpp

#include "SimdArray.hpp"

namespace t2chrr { // t2chrr namespace

/// @brief Computes contracted [F|X|D]  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fd The index of integral in contracted integrals buffer.
/// @param idx_fp The index of integral in contracted integrals buffer.
/// @param idx_gp The index of integral in contracted integrals buffer.
/// @param factors The contracted factors buffer.
auto
comp_hrr_fd(CSimdArray<double>& cbuffer, 
            const size_t idx_fd,
            const size_t idx_fp,
            const size_t idx_gp,
            const CSimdArray<double>& factors) -> void;
} // t2chrr namespace

#endif /* T2CHrrABRecFD_hpp */

