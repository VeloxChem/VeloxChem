#ifndef ElectricDipoleMomentumPrimRecHH
#define ElectricDipoleMomentumPrimRecHH

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [H|r|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_hh The index of integral in primitive integrals buffer.
/// @param idx_dip_fh The index of integral in primitive integrals buffer.
/// @param idx_dip_gg The index of integral in primitive integrals buffer.
/// @param idx_ovl_gh The index of integral in primitive integrals buffer.
/// @param idx_dip_gh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_hh(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_hh,
                                           const size_t              idx_dip_fh,
                                           const size_t              idx_dip_gg,
                                           const size_t              idx_ovl_gh,
                                           const size_t              idx_dip_gh,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecHH */
