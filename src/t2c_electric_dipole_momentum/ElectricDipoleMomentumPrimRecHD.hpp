#ifndef ElectricDipoleMomentumPrimRecHD
#define ElectricDipoleMomentumPrimRecHD

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [H|r|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_hd The index of integral in primitive integrals buffer.
/// @param idx_dip_fd The index of integral in primitive integrals buffer.
/// @param idx_dip_gp The index of integral in primitive integrals buffer.
/// @param idx_ovl_gd The index of integral in primitive integrals buffer.
/// @param idx_dip_gd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_hd(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_hd,
                                      const size_t idx_dip_fd,
                                      const size_t idx_dip_gp,
                                      const size_t idx_ovl_gd,
                                      const size_t idx_dip_gd,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecHD */
