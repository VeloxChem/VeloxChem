#ifndef ThreeCenterOverlapPrimRecFF
#define ThreeCenterOverlapPrimRecFF

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [F|G(r)|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ff The index of integral in primitive integrals buffer.
/// @param idx_pf The index of integral in primitive integrals buffer.
/// @param idx_dd The index of integral in primitive integrals buffer.
/// @param idx_df The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_ff(CSimdArray<double>& pbuffer, 
                     const size_t idx_ff,
                     const size_t idx_pf,
                     const size_t idx_dd,
                     const size_t idx_df,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecFF */
