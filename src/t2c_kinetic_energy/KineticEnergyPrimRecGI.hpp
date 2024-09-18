#ifndef KineticEnergyPrimRecGI
#define KineticEnergyPrimRecGI

#include "SimdArray.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes primitive [G|T|I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_gi The index of integral in primitive integrals buffer.
/// @param idx_ovl_di The index of integral in primitive integrals buffer.
/// @param idx_kin_di The index of integral in primitive integrals buffer.
/// @param idx_kin_fh The index of integral in primitive integrals buffer.
/// @param idx_kin_fi The index of integral in primitive integrals buffer.
/// @param idx_ovl_gi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto comp_prim_kinetic_energy_gi(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_kin_gi,
                                 const size_t              idx_ovl_di,
                                 const size_t              idx_kin_di,
                                 const size_t              idx_kin_fh,
                                 const size_t              idx_kin_fi,
                                 const size_t              idx_ovl_gi,
                                 const CSimdArray<double>& factors,
                                 const size_t              idx_rpa,
                                 const double              a_exp) -> void;
}  // namespace kinrec

#endif /* KineticEnergyPrimRecGI */