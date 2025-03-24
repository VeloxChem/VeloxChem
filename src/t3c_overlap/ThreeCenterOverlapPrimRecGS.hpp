#ifndef ThreeCenterOverlapPrimRecGS
#define ThreeCenterOverlapPrimRecGS

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [G|G(r)|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gs The index of integral in primitive integrals buffer.
/// @param idx_ds The index of integral in primitive integrals buffer.
/// @param idx_fs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gs(CSimdArray<double>& pbuffer, 
                     const size_t idx_gs,
                     const size_t idx_ds,
                     const size_t idx_fs,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecGS */
