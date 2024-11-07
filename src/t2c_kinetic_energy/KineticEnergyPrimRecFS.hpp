#ifndef KineticEnergyPrimRecFS
#define KineticEnergyPrimRecFS

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [F|T|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_fs The index of integral in primitive integrals buffer.
/// @param idx_ovl_ps The index of integral in primitive integrals buffer.
/// @param idx_kin_ps The index of integral in primitive integrals buffer.
/// @param idx_kin_ds The index of integral in primitive integrals buffer.
/// @param idx_ovl_fs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_fs(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_fs,
                                 const size_t              idx_ovl_ps,
                                 const size_t              idx_kin_ps,
                                 const size_t              idx_kin_ds,
                                 const size_t              idx_ovl_fs,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecFS */
