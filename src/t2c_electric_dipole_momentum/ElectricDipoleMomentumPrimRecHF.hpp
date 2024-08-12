#ifndef ElectricDipoleMomentumPrimRecHF
#define ElectricDipoleMomentumPrimRecHF

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [H|r|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_hf The index of integral in primitive integrals buffer.
/// @param idx_dip_ff The index of integral in primitive integrals buffer.
/// @param idx_dip_gd The index of integral in primitive integrals buffer.
/// @param idx_ovl_gf The index of integral in primitive integrals buffer.
/// @param idx_dip_gf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_hf(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_hf,
                                           const size_t              idx_dip_ff,
                                           const size_t              idx_dip_gd,
                                           const size_t              idx_ovl_gf,
                                           const size_t              idx_dip_gf,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecHF */
