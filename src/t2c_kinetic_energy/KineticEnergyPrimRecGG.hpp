#ifndef KineticEnergyPrimRecGG
#define KineticEnergyPrimRecGG

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [G|T|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_gg The index of integral in primitive integrals buffer.
/// @param idx_ovl_dg The index of integral in primitive integrals buffer.
/// @param idx_kin_dg The index of integral in primitive integrals buffer.
/// @param idx_kin_ff The index of integral in primitive integrals buffer.
/// @param idx_kin_fg The index of integral in primitive integrals buffer.
/// @param idx_ovl_gg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_gg(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_gg,
                                 const size_t              idx_ovl_dg,
                                 const size_t              idx_kin_dg,
                                 const size_t              idx_kin_ff,
                                 const size_t              idx_kin_fg,
                                 const size_t              idx_ovl_gg,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecGG */
