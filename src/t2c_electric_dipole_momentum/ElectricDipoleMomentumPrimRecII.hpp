#ifndef ElectricDipoleMomentumPrimRecII
#define ElectricDipoleMomentumPrimRecII

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [I|r|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_ii The index of integral in primitive integrals buffer.
/// @param idx_dip_gi The index of integral in primitive integrals buffer.
/// @param idx_dip_hh The index of integral in primitive integrals buffer.
/// @param idx_ovl_hi The index of integral in primitive integrals buffer.
/// @param idx_dip_hi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_electric_dipole_momentum_ii(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_ii,
                                           const size_t              idx_dip_gi,
                                           const size_t              idx_dip_hh,
                                           const size_t              idx_ovl_hi,
                                           const size_t              idx_dip_hi,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpa,
                                           const double              a_exp) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecII */
