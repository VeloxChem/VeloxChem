#ifndef ThreeCenterOverlapPrimRecGF
#define ThreeCenterOverlapPrimRecGF

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [G|G(r)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_gf The index of integral in primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param idx_fd The index of integral in primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_gf(CSimdArray<double>& pbuffer, 
                     const size_t idx_gf,
                     const size_t idx_df,
                     const size_t idx_fd,
                     const size_t idx_ff,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecGF */
