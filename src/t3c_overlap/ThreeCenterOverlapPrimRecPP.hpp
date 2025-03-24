#ifndef ThreeCenterOverlapPrimRecPP
#define ThreeCenterOverlapPrimRecPP

#include "SimdArray.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes primitive [P|G(r)|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_pp The index of integral in primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param idx_sp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rga The vector of distances R(GA) = G - A.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The primitive basis function exponent on center C.
auto
comp_prim_overlap_pp(CSimdArray<double>& pbuffer, 
                     const size_t idx_pp,
                     const size_t idx_ss,
                     const size_t idx_sp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void;
} // t3ovlrec namespace

#endif /* ThreeCenterOverlapPrimRecPP */
