#ifndef ElectricDipoleMomentumPrimRecGH
#define ElectricDipoleMomentumPrimRecGH

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [G|r|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_gh The index of integral in primitive integrals buffer.
/// @param idx_dip_dh The index of integral in primitive integrals buffer.
/// @param idx_dip_fg The index of integral in primitive integrals buffer.
/// @param idx_ovl_fh The index of integral in primitive integrals buffer.
/// @param idx_dip_fh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_gh(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_gh,
                                      const size_t idx_dip_dh,
                                      const size_t idx_dip_fg,
                                      const size_t idx_ovl_fh,
                                      const size_t idx_dip_fh,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecGH */
