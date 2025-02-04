#ifndef TwoCenterElectronRepulsionPrimRecFI
#define TwoCenterElectronRepulsionPrimRecFI

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [F|1/|r-r'||I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_pi The index of integral in primitive integrals buffer.
/// @param idx_eri_1_pi The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_di The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_fi(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fi,
                                const size_t idx_eri_0_pi,
                                const size_t idx_eri_1_pi,
                                const size_t idx_eri_1_dh,
                                const size_t idx_eri_1_di,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecFI */
