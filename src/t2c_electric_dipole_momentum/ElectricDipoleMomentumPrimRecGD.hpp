#ifndef ElectricDipoleMomentumPrimRecGD
#define ElectricDipoleMomentumPrimRecGD

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [G|r|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_gd The index of integral in primitive integrals buffer.
/// @param idx_dip_dd The index of integral in primitive integrals buffer.
/// @param idx_dip_fp The index of integral in primitive integrals buffer.
/// @param idx_ovl_fd The index of integral in primitive integrals buffer.
/// @param idx_dip_fd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_gd(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_gd,
                                           const size_t              idx_dip_dd,
                                           const size_t              idx_dip_fp,
                                           const size_t              idx_ovl_fd,
                                           const size_t              idx_dip_fd,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecGD */
