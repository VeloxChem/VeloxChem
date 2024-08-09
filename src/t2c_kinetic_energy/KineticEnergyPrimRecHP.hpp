#ifndef KineticEnergyPrimRecHP
#define KineticEnergyPrimRecHP

#include "SimdArray.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes primitive [H|T|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_hp The index of integral in primitive integrals buffer.
/// @param idx_ovl_fp The index of integral in primitive integrals buffer.
/// @param idx_kin_fp The index of integral in primitive integrals buffer.
/// @param idx_kin_gs The index of integral in primitive integrals buffer.
/// @param idx_kin_gp The index of integral in primitive integrals buffer.
/// @param idx_ovl_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_kinetic_energy_hp(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_hp,
                            const size_t idx_ovl_fp,
                            const size_t idx_kin_fp,
                            const size_t idx_kin_gs,
                            const size_t idx_kin_gp,
                            const size_t idx_ovl_hp,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void;
} // kinrec namespace

#endif /* KineticEnergyPrimRecHP */
