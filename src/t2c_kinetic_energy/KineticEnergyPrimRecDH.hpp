#ifndef KineticEnergyPrimRecDH
#define KineticEnergyPrimRecDH

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [D|T|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_dh The index of integral in primitive integrals buffer.
/// @param idx_ovl_sh The index of integral in primitive integrals buffer.
/// @param idx_kin_sh The index of integral in primitive integrals buffer.
/// @param idx_kin_pg The index of integral in primitive integrals buffer.
/// @param idx_kin_ph The index of integral in primitive integrals buffer.
/// @param idx_ovl_dh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_dh(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_dh,
                                 const size_t              idx_ovl_sh,
                                 const size_t              idx_kin_sh,
                                 const size_t              idx_kin_pg,
                                 const size_t              idx_kin_ph,
                                 const size_t              idx_ovl_dh,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecDH */
