#ifndef ThreeCenterOverlapPrimRecSS_hpp
#define ThreeCenterOverlapPrimRecSS_hpp

#include "SimdArray.hpp"

namespace t3ovlrec {  // t3ovlrec namespace

/// @brief Computes primitive [S|A|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_t3ovl_ss The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss  The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpc The vector of distances R(PC) = P - C.
/// @param a_exp The primitive basis function exponent on center A.
/// @param g_exp The external basis function exponent on center C.
/// @param g_norm The external basis function normalization factor on center C.
auto comp_prim_overlap_ss(CSimdArray<double>&       pbuffer,
                          const size_t              idx_t3ovl_ss,
                          const size_t              idx_ovl_ss,
                          CSimdArray<double>&       factors,
                          const size_t              idx_rpc, 
                          const double              a_exp,
                          const double              g_exp,
                          const double              g_norm) -> void;
}  // namespace npotrec

#endif /* ThreeCenterOverlapPrimRecSS_hpp */
