#ifndef KineticEnergyPrimRecIH
#define KineticEnergyPrimRecIH

#include "SimdArray.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes primitive [I|T|H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_ih The index of integral in primitive integrals buffer.
/// @param idx_ovl_gh The index of integral in primitive integrals buffer.
/// @param idx_kin_gh The index of integral in primitive integrals buffer.
/// @param idx_kin_hg The index of integral in primitive integrals buffer.
/// @param idx_kin_hh The index of integral in primitive integrals buffer.
/// @param idx_ovl_ih The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_kinetic_energy_ih(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_ih,
                            const size_t idx_ovl_gh,
                            const size_t idx_kin_gh,
                            const size_t idx_kin_hg,
                            const size_t idx_kin_hh,
                            const size_t idx_ovl_ih,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void;
} // kinrec namespace

#endif /* KineticEnergyPrimRecIH */
