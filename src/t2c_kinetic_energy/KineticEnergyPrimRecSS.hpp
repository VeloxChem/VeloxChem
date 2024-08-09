#ifndef KineticEnergyPrimRecSS
#define KineticEnergyPrimRecSS

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [S|T|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_ss The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_ss(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_ss,
                                 const size_t              idx_ovl_ss,
                                 const CSimdArray<double>& factors,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecSS */
