#ifndef ElectricDipoleMomentumPrimRecPS
#define ElectricDipoleMomentumPrimRecPS

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [P|r|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_ps The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param idx_dip_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_ps(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_ps,
                                           const size_t              idx_ovl_ss,
                                           const size_t              idx_dip_ss,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecPS */
