#ifndef ElectricDipoleMomentumPrimRecFS
#define ElectricDipoleMomentumPrimRecFS

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [F|r|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_fs The index of integral in primitive integrals buffer.
/// @param idx_dip_ps The index of integral in primitive integrals buffer.
/// @param idx_ovl_ds The index of integral in primitive integrals buffer.
/// @param idx_dip_ds The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_fs(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_fs,
                                      const size_t idx_dip_ps,
                                      const size_t idx_ovl_ds,
                                      const size_t idx_dip_ds,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecFS */
