#ifndef ElectricDipoleMomentumPrimRecPF
#define ElectricDipoleMomentumPrimRecPF

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [P|r|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_pf The index of integral in primitive integrals buffer.
/// @param idx_dip_sd The index of integral in primitive integrals buffer.
/// @param idx_ovl_sf The index of integral in primitive integrals buffer.
/// @param idx_dip_sf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_pf(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_pf,
                                           const size_t              idx_dip_sd,
                                           const size_t              idx_ovl_sf,
                                           const size_t              idx_dip_sf,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecPF */
