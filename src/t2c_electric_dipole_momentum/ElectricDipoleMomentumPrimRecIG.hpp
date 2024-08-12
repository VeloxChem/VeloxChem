#ifndef ElectricDipoleMomentumPrimRecIG
#define ElectricDipoleMomentumPrimRecIG

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [I|r|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_ig The index of integral in primitive integrals buffer.
/// @param idx_dip_gg The index of integral in primitive integrals buffer.
/// @param idx_dip_hf The index of integral in primitive integrals buffer.
/// @param idx_ovl_hg The index of integral in primitive integrals buffer.
/// @param idx_dip_hg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_ig(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_ig,
                                           const size_t              idx_dip_gg,
                                           const size_t              idx_dip_hf,
                                           const size_t              idx_ovl_hg,
                                           const size_t              idx_dip_hg,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecIG */
