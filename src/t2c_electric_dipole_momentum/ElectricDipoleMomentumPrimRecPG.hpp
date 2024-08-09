#ifndef ElectricDipoleMomentumPrimRecPG
#define ElectricDipoleMomentumPrimRecPG

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [P|r|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_pg The index of integral in primitive integrals buffer.
/// @param idx_dip_sf The index of integral in primitive integrals buffer.
/// @param idx_ovl_sg The index of integral in primitive integrals buffer.
/// @param idx_dip_sg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_pg(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_pg,
                                      const size_t idx_dip_sf,
                                      const size_t idx_ovl_sg,
                                      const size_t idx_dip_sg,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecPG */
