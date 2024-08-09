#ifndef ElectricDipoleMomentumPrimRecGS
#define ElectricDipoleMomentumPrimRecGS

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [G|r|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_gs The index of integral in primitive integrals buffer.
/// @param idx_dip_ds The index of integral in primitive integrals buffer.
/// @param idx_ovl_fs The index of integral in primitive integrals buffer.
/// @param idx_dip_fs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_gs(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_gs,
                                      const size_t idx_dip_ds,
                                      const size_t idx_ovl_fs,
                                      const size_t idx_dip_fs,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecGS */
