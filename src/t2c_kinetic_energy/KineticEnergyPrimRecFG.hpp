#ifndef KineticEnergyPrimRecFG
#define KineticEnergyPrimRecFG

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [F|T|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_fg The index of integral in primitive integrals buffer.
/// @param idx_ovl_pg The index of integral in primitive integrals buffer.
/// @param idx_kin_pg The index of integral in primitive integrals buffer.
/// @param idx_kin_df The index of integral in primitive integrals buffer.
/// @param idx_kin_dg The index of integral in primitive integrals buffer.
/// @param idx_ovl_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_fg(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_fg,
                                 const size_t              idx_ovl_pg,
                                 const size_t              idx_kin_pg,
                                 const size_t              idx_kin_df,
                                 const size_t              idx_kin_dg,
                                 const size_t              idx_ovl_fg,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecFG */
