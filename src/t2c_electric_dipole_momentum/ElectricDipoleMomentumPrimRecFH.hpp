#ifndef ElectricDipoleMomentumPrimRecFH
#define ElectricDipoleMomentumPrimRecFH

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [F|r|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_fh The index of integral in primitive integrals buffer.
/// @param idx_dip_ph The index of integral in primitive integrals buffer.
/// @param idx_dip_dg The index of integral in primitive integrals buffer.
/// @param idx_ovl_dh The index of integral in primitive integrals buffer.
/// @param idx_dip_dh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_fh(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_fh,
                                      const size_t idx_dip_ph,
                                      const size_t idx_dip_dg,
                                      const size_t idx_ovl_dh,
                                      const size_t idx_dip_dh,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecFH */
