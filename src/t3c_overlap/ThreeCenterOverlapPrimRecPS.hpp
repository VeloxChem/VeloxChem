#ifndef ThreeCenterOverlapPrimRecPS
#define ThreeCenterOverlapPrimRecPS

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [P|G(r)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ps The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
auto
comp_prim_overlap_ps(CSimdArray<double>& pbuffer, 
                     const size_t idx_ps,
                     const size_t idx_ss,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecPS */
