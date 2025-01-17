#ifndef TwoCenterElectronRepulsionPrimRecDS
#define TwoCenterElectronRepulsionPrimRecDS

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [D|1/|r-r'||S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ds The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ss The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ss The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ps The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_ds(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ds,
                                const size_t idx_eri_0_ss,
                                const size_t idx_eri_1_ss,
                                const size_t idx_eri_1_ps,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecDS */
