#ifndef KineticEnergyPrimRecDP
#define KineticEnergyPrimRecDP

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [D|T|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_dp The index of integral in primitive integrals buffer.
/// @param idx_ovl_sp The index of integral in primitive integrals buffer.
/// @param idx_kin_sp The index of integral in primitive integrals buffer.
/// @param idx_kin_ps The index of integral in primitive integrals buffer.
/// @param idx_kin_pp The index of integral in primitive integrals buffer.
/// @param idx_ovl_dp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_dp(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_dp,
                                 const size_t              idx_ovl_sp,
                                 const size_t              idx_kin_sp,
                                 const size_t              idx_kin_ps,
                                 const size_t              idx_kin_pp,
                                 const size_t              idx_ovl_dp,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecDP */
