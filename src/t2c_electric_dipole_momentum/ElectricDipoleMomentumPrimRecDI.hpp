#ifndef ElectricDipoleMomentumPrimRecDI
#define ElectricDipoleMomentumPrimRecDI

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [D|r|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_di The index of integral in primitive integrals buffer.
/// @param idx_dip_si The index of integral in primitive integrals buffer.
/// @param idx_dip_ph The index of integral in primitive integrals buffer.
/// @param idx_ovl_pi The index of integral in primitive integrals buffer.
/// @param idx_dip_pi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_di(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_di,
                                      const size_t idx_dip_si,
                                      const size_t idx_dip_ph,
                                      const size_t idx_ovl_pi,
                                      const size_t idx_dip_pi,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecDI */
