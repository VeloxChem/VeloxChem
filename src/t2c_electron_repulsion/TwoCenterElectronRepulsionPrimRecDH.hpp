#ifndef TwoCenterElectronRepulsionPrimRecDH
#define TwoCenterElectronRepulsionPrimRecDH

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [D|1/|r-r'||H]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_pg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ph The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_dh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_dh,
                                const size_t idx_eri_0_sh,
                                const size_t idx_eri_1_sh,
                                const size_t idx_eri_1_pg,
                                const size_t idx_eri_1_ph,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecDH */
