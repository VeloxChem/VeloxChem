#ifndef TwoCenterElectronRepulsionPrimRecDI
#define TwoCenterElectronRepulsionPrimRecDI

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [D|1/|r-r'||I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_di The index of integral in primitive integrals buffer.
/// @param idx_eri_0_si The index of integral in primitive integrals buffer.
/// @param idx_eri_1_si The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ph The index of integral in primitive integrals buffer.
/// @param idx_eri_1_pi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_di(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_di,
                                const size_t idx_eri_0_si,
                                const size_t idx_eri_1_si,
                                const size_t idx_eri_1_ph,
                                const size_t idx_eri_1_pi,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecDI */
