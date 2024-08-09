#ifndef ElectricDipoleMomentumPrimRecSP
#define ElectricDipoleMomentumPrimRecSP

#include "SimdArray.hpp"

namespace diprec { // diprec namespace

/// @brief Computes primitive [S|r|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_sp The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param idx_dip_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electric_dipole_momentum_sp(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_sp,
                                      const size_t idx_ovl_ss,
                                      const size_t idx_dip_ss,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpb,
                                      const double a_exp) -> void;
} // diprec namespace

#endif /* ElectricDipoleMomentumPrimRecSP */
