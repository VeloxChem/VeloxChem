#ifndef TwoCenterElectronRepulsionPrimRecFH
#define TwoCenterElectronRepulsionPrimRecFH

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [F|1/|r-r'||H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ph The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ph The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dh The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_fh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fh,
                                const size_t idx_eri_0_ph,
                                const size_t idx_eri_1_ph,
                                const size_t idx_eri_1_dg,
                                const size_t idx_eri_1_dh,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecFH */
