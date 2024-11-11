#ifndef KineticEnergyPrimRecIF
#define KineticEnergyPrimRecIF

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [I|T|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_if The index of integral in primitive integrals buffer.
/// @param idx_ovl_gf The index of integral in primitive integrals buffer.
/// @param idx_kin_gf The index of integral in primitive integrals buffer.
/// @param idx_kin_hd The index of integral in primitive integrals buffer.
/// @param idx_kin_hf The index of integral in primitive integrals buffer.
/// @param idx_ovl_if The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_if(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_if,
                                 const size_t              idx_ovl_gf,
                                 const size_t              idx_kin_gf,
                                 const size_t              idx_kin_hd,
                                 const size_t              idx_kin_hf,
                                 const size_t              idx_ovl_if,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecIF */
