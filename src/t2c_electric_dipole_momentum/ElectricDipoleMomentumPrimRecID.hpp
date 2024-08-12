#ifndef ElectricDipoleMomentumPrimRecID
#define ElectricDipoleMomentumPrimRecID

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [I|r|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_id The index of integral in primitive integrals buffer.
/// @param idx_dip_gd The index of integral in primitive integrals buffer.
/// @param idx_dip_hp The index of integral in primitive integrals buffer.
/// @param idx_ovl_hd The index of integral in primitive integrals buffer.
/// @param idx_dip_hd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_id(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_id,
                                      const size_t idx_dip_gd,
                                      const size_t idx_dip_hp,
                                      const size_t idx_ovl_hd,
                                      const size_t idx_dip_hd,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecID */
