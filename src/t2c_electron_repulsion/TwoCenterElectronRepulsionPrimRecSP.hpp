#ifndef TwoCenterElectronRepulsionPrimRecSP
#define TwoCenterElectronRepulsionPrimRecSP

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [S|1/|r-r'||P]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
auto
comp_prim_electron_repulsion_sp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sp,
                                const size_t idx_eri_1_ss,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecSP */
