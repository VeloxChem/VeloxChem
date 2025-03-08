#ifndef TwoCenterElectronRepulsionPrimRecPS
#define TwoCenterElectronRepulsionPrimRecPS

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [P|1/|r-r'||S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ps The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
auto
comp_prim_electron_repulsion_ps(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ps,
                                const size_t idx_eri_1_ss,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecPS */
