#ifndef KineticEnergyPrimRecHS
#define KineticEnergyPrimRecHS

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [H|T|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_hs The index of integral in primitive integrals buffer.
/// @param idx_ovl_fs The index of integral in primitive integrals buffer.
/// @param idx_kin_fs The index of integral in primitive integrals buffer.
/// @param idx_kin_gs The index of integral in primitive integrals buffer.
/// @param idx_ovl_hs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_hs(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_hs,
                                 const size_t              idx_ovl_fs,
                                 const size_t              idx_kin_fs,
                                 const size_t              idx_kin_gs,
                                 const size_t              idx_ovl_hs,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecHS */
