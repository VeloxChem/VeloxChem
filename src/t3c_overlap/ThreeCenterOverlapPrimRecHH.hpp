#ifndef ThreeCenterOverlapPrimRecHH
#define ThreeCenterOverlapPrimRecHH

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [H|G(r)|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_hh The index of integral in primitive integrals buffer.
/// @param idx_fh The index of integral in primitive integrals buffer.
/// @param idx_gg The index of integral in primitive integrals buffer.
/// @param idx_gh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_hh(CSimdArray<double>& pbuffer, 
                     const size_t idx_hh,
                     const size_t idx_fh,
                     const size_t idx_gg,
                     const size_t idx_gh,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecHH */
