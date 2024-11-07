#ifndef OverlapPrimRecSS
#define OverlapPrimRecSS

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [S|S]  integrals for set of data buffers.
/// @param pbuffer  The primitive integrals buffer.
/// @param idx_ovl_ss The index of integrals in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
/// @param a_norm The primitive basis function normalization factor on center A.
auto comp_prim_overlap_ss(CSimdArray<double>& pbuffer, const size_t idx_ovl_ss, CSimdArray<double>& factors, const double a_exp, const double a_norm)
    -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecSS */
