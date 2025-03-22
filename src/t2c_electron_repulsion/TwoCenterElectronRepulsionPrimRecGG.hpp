#ifndef TwoCenterElectronRepulsionPrimRecGG
#define TwoCenterElectronRepulsionPrimRecGG

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [G|1/|r-r'||G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ff The index of integral in primitive integrals buffer.
/// @param idx_eri_1_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_gg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gg,
                                const size_t idx_eri_0_dg,
                                const size_t idx_eri_1_dg,
                                const size_t idx_eri_1_ff,
                                const size_t idx_eri_1_fg,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecGG */
