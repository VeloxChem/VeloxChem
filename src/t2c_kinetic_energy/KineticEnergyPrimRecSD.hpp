#ifndef KineticEnergyPrimRecSD
#define KineticEnergyPrimRecSD

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [S|T|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_sd The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param idx_kin_ss The index of integral in primitive integrals buffer.
/// @param idx_kin_sp The index of integral in primitive integrals buffer.
/// @param idx_ovl_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_sd(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_sd,
                                 const size_t              idx_ovl_ss,
                                 const size_t              idx_kin_ss,
                                 const size_t              idx_kin_sp,
                                 const size_t              idx_ovl_sd,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpb,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecSD */
