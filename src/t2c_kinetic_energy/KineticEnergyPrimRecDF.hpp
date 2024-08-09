#ifndef KineticEnergyPrimRecDF
#define KineticEnergyPrimRecDF

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [D|T|F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_df The index of integral in primitive integrals buffer.
/// @param idx_ovl_sf The index of integral in primitive integrals buffer.
/// @param idx_kin_sf The index of integral in primitive integrals buffer.
/// @param idx_kin_pd The index of integral in primitive integrals buffer.
/// @param idx_kin_pf The index of integral in primitive integrals buffer.
/// @param idx_ovl_df The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_df(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_df,
                                 const size_t              idx_ovl_sf,
                                 const size_t              idx_kin_sf,
                                 const size_t              idx_kin_pd,
                                 const size_t              idx_kin_pf,
                                 const size_t              idx_ovl_df,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecDF */
