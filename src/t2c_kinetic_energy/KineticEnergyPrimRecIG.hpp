#ifndef KineticEnergyPrimRecIG
#define KineticEnergyPrimRecIG

#include "SimdArray.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes primitive [I|T|G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_kin_ig The index of integral in primitive integrals buffer.
/// @param idx_ovl_gg The index of integral in primitive integrals buffer.
/// @param idx_kin_gg The index of integral in primitive integrals buffer.
/// @param idx_kin_hf The index of integral in primitive integrals buffer.
/// @param idx_kin_hg The index of integral in primitive integrals buffer.
/// @param idx_ovl_ig The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_kinetic_energy_ig(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_ig,
                            const size_t idx_ovl_gg,
                            const size_t idx_kin_gg,
                            const size_t idx_kin_hf,
                            const size_t idx_kin_hg,
                            const size_t idx_ovl_ig,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpa,
                            const double a_exp) -> void;
} // kinrec namespace

#endif /* KineticEnergyPrimRecIG */
