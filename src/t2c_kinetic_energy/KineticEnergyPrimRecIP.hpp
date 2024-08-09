#ifndef KineticEnergyPrimRecIP
#define KineticEnergyPrimRecIP

#include "SimdArray.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes primitive [I|T|P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_ip The index of integral in primitive integrals buffer.
/// @param idx_ovl_gp The index of integral in primitive integrals buffer.
/// @param idx_kin_gp The index of integral in primitive integrals buffer.
/// @param idx_kin_hs The index of integral in primitive integrals buffer.
/// @param idx_kin_hp The index of integral in primitive integrals buffer.
/// @param idx_ovl_ip The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_kinetic_energy_ip(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_ip,
                            const size_t idx_ovl_gp,
                            const size_t idx_kin_gp,
                            const size_t idx_kin_hs,
                            const size_t idx_kin_hp,
                            const size_t idx_ovl_ip,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void;
} // kinrec namespace

#endif /* KineticEnergyPrimRecIP */
