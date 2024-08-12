#ifndef ElectricDipoleMomentumPrimRecIP
#define ElectricDipoleMomentumPrimRecIP

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [I|r|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_ip The index of integral in primitive integrals buffer.
/// @param idx_dip_gp The index of integral in primitive integrals buffer.
/// @param idx_dip_hs The index of integral in primitive integrals buffer.
/// @param idx_ovl_hp The index of integral in primitive integrals buffer.
/// @param idx_dip_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_ip(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ip,
                                      const size_t idx_dip_gp,
                                      const size_t idx_dip_hs,
                                      const size_t idx_ovl_hp,
                                      const size_t idx_dip_hp,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecIP */
