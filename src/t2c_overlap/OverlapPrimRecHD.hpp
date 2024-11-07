#ifndef OverlapPrimRecHD
#define OverlapPrimRecHD

#include "SimdArray.hpp"

namespace ovlrec {  // ovlrec namespace

/// @brief Computes primitive [H|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_hd The index of integral in primitive integrals buffer.
/// @param idx_ovl_fd The index of integral in primitive integrals buffer.
/// @param idx_ovl_gp The index of integral in primitive integrals buffer.
/// @param idx_ovl_gd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_overlap_hd(CSimdArray<double>&       pbuffer,
                          const size_t              idx_ovl_hd,
                          const size_t              idx_ovl_fd,
                          const size_t              idx_ovl_gp,
                          const size_t              idx_ovl_gd,
                          const CSimdArray<double>& factors,
                          const size_t              idx_rpa,
                          const double              a_exp) -> void;
}  // namespace ovlrec

#endif /* OverlapPrimRecHD */
