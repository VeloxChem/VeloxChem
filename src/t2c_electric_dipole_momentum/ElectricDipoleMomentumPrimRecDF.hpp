#ifndef ElectricDipoleMomentumPrimRecDF
#define ElectricDipoleMomentumPrimRecDF

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [D|r|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_df The index of integral in primitive integrals buffer.
/// @param idx_dip_sf The index of integral in primitive integrals buffer.
/// @param idx_dip_pd The index of integral in primitive integrals buffer.
/// @param idx_ovl_pf The index of integral in primitive integrals buffer.
/// @param idx_dip_pf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_df(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_df,
                                           const size_t              idx_dip_sf,
                                           const size_t              idx_dip_pd,
                                           const size_t              idx_ovl_pf,
                                           const size_t              idx_dip_pf,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecDF */
