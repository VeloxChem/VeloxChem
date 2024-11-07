#ifndef KineticEnergyPrimRecPF
#define KineticEnergyPrimRecPF

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [P|T|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_pf The index of integral in primitive integrals buffer.
/// @param idx_kin_sd The index of integral in primitive integrals buffer.
/// @param idx_kin_sf The index of integral in primitive integrals buffer.
/// @param idx_ovl_pf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_pf(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_pf,
                                 const size_t              idx_kin_sd,
                                 const size_t              idx_kin_sf,
                                 const size_t              idx_ovl_pf,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecPF */
