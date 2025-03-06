#ifndef OverlapPrimRecSI
#define OverlapPrimRecSI

#include "SimdArray.hpp"

namespace ovlrec { // ovlrec namespace

/// @brief Computes primitive [S|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ovl_si The index of integral in primitive integrals buffer.
/// @param idx_ovl_sg The index of integral in primitive integrals buffer.
/// @param idx_ovl_sh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_overlap_si(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_si,
                     const size_t idx_ovl_sg,
                     const size_t idx_ovl_sh,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpb,
                     const double a_exp) -> void;
} // ovlrec namespace

#endif /* OverlapPrimRecSI */
