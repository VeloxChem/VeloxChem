#ifndef KineticEnergyPrimRecHI
#define KineticEnergyPrimRecHI

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [H|T|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_hi The index of integral in primitive integrals buffer.
/// @param idx_ovl_fi The index of integral in primitive integrals buffer.
/// @param idx_kin_fi The index of integral in primitive integrals buffer.
/// @param idx_kin_gh The index of integral in primitive integrals buffer.
/// @param idx_kin_gi The index of integral in primitive integrals buffer.
/// @param idx_ovl_hi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_hi(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_hi,
                                 const size_t              idx_ovl_fi,
                                 const size_t              idx_kin_fi,
                                 const size_t              idx_kin_gh,
                                 const size_t              idx_kin_gi,
                                 const size_t              idx_ovl_hi,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecHI */
