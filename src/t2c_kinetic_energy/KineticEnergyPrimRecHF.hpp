#ifndef KineticEnergyPrimRecHF
#define KineticEnergyPrimRecHF

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [H|T|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_hf The index of integral in primitive integrals buffer.
/// @param idx_ovl_ff The index of integral in primitive integrals buffer.
/// @param idx_kin_ff The index of integral in primitive integrals buffer.
/// @param idx_kin_gd The index of integral in primitive integrals buffer.
/// @param idx_kin_gf The index of integral in primitive integrals buffer.
/// @param idx_ovl_hf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_hf(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_hf,
                                 const size_t              idx_ovl_ff,
                                 const size_t              idx_kin_ff,
                                 const size_t              idx_kin_gd,
                                 const size_t              idx_kin_gf,
                                 const size_t              idx_ovl_hf,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecHF */
