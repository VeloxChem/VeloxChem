#ifndef ElectricDipoleMomentumPrimRecFI
#define ElectricDipoleMomentumPrimRecFI

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [F|r|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_fi The index of integral in primitive integrals buffer.
/// @param idx_dip_pi The index of integral in primitive integrals buffer.
/// @param idx_dip_dh The index of integral in primitive integrals buffer.
/// @param idx_ovl_di The index of integral in primitive integrals buffer.
/// @param idx_dip_di The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_fi(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_fi,
                                      const size_t idx_dip_pi,
                                      const size_t idx_dip_dh,
                                      const size_t idx_ovl_di,
                                      const size_t idx_dip_di,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecFI */
