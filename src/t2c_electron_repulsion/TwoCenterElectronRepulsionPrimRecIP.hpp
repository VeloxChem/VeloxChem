#ifndef TwoCenterElectronRepulsionPrimRecIP
#define TwoCenterElectronRepulsionPrimRecIP

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [I|1/|r-r'||P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ip The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_gp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hs The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_ip(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ip,
                                const size_t idx_eri_0_gp,
                                const size_t idx_eri_1_gp,
                                const size_t idx_eri_1_hs,
                                const size_t idx_eri_1_hp,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecIP */
