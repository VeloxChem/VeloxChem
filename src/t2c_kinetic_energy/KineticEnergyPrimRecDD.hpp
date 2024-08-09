#ifndef KineticEnergyPrimRecDD
#define KineticEnergyPrimRecDD

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [D|T|D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_dd The index of integral in primitive integrals buffer.
/// @param idx_ovl_sd The index of integral in primitive integrals buffer.
/// @param idx_kin_sd The index of integral in primitive integrals buffer.
/// @param idx_kin_pp The index of integral in primitive integrals buffer.
/// @param idx_kin_pd The index of integral in primitive integrals buffer.
/// @param idx_ovl_dd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_dd(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_dd,
                                 const size_t              idx_ovl_sd,
                                 const size_t              idx_kin_sd,
                                 const size_t              idx_kin_pp,
                                 const size_t              idx_kin_pd,
                                 const size_t              idx_ovl_dd,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecDD */
