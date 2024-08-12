#ifndef ElectricDipoleMomentumPrimRecHG
#define ElectricDipoleMomentumPrimRecHG

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [H|r|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_hg The index of integral in primitive integrals buffer.
/// @param idx_dip_fg The index of integral in primitive integrals buffer.
/// @param idx_dip_gf The index of integral in primitive integrals buffer.
/// @param idx_ovl_gg The index of integral in primitive integrals buffer.
/// @param idx_dip_gg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_hg(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_hg,
                                           const size_t              idx_dip_fg,
                                           const size_t              idx_dip_gf,
                                           const size_t              idx_ovl_gg,
                                           const size_t              idx_dip_gg,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecHG */
