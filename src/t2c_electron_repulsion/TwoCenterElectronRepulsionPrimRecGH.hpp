#ifndef TwoCenterElectronRepulsionPrimRecGH
#define TwoCenterElectronRepulsionPrimRecGH

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [G|1/|r-r'||H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_fg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_fh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_gh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gh,
                                const size_t idx_eri_0_dh,
                                const size_t idx_eri_1_dh,
                                const size_t idx_eri_1_fg,
                                const size_t idx_eri_1_fh,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecGH */
