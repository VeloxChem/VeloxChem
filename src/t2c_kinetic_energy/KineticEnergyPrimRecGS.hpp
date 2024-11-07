#ifndef KineticEnergyPrimRecGS
#define KineticEnergyPrimRecGS

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [G|T|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_gs The index of integral in primitive integrals buffer.
/// @param idx_ovl_ds The index of integral in primitive integrals buffer.
/// @param idx_kin_ds The index of integral in primitive integrals buffer.
/// @param idx_kin_fs The index of integral in primitive integrals buffer.
/// @param idx_ovl_gs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_gs(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_gs,
                                 const size_t              idx_ovl_ds,
                                 const size_t              idx_kin_ds,
                                 const size_t              idx_kin_fs,
                                 const size_t              idx_ovl_gs,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecGS */
