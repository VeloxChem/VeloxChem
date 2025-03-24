#ifndef ThreeCenterOverlapPrimRecSP
#define ThreeCenterOverlapPrimRecSP

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [S|G(r)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rgb The vector of distances R(GB) = G - B.
auto
comp_prim_overlap_sp(CSimdArray<double>& pbuffer, 
                     const size_t idx_sp,
                     const size_t idx_ss,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecSP */
