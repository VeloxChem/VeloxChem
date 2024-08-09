#ifndef KineticEnergyPrimRecHG
#define KineticEnergyPrimRecHG

#include "SimdArray.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes primitive [H|T|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_hg The index of integral in primitive integrals buffer.
/// @param idx_ovl_fg The index of integral in primitive integrals buffer.
/// @param idx_kin_fg The index of integral in primitive integrals buffer.
/// @param idx_kin_gf The index of integral in primitive integrals buffer.
/// @param idx_kin_gg The index of integral in primitive integrals buffer.
/// @param idx_ovl_hg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_kinetic_energy_hg(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_hg,
                            const size_t idx_ovl_fg,
                            const size_t idx_kin_fg,
                            const size_t idx_kin_gf,
                            const size_t idx_kin_gg,
                            const size_t idx_ovl_hg,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void;
} // kinrec namespace

#endif /* KineticEnergyPrimRecHG */
