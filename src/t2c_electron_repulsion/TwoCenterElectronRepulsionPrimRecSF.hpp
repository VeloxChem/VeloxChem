#ifndef TwoCenterElectronRepulsionPrimRecSF
#define TwoCenterElectronRepulsionPrimRecSF

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [S|1/|r-r'||F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpb The vector of distances R(PB) = P - B.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_sf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sf,
                                const size_t idx_eri_0_sp,
                                const size_t idx_eri_1_sp,
                                const size_t idx_eri_1_sd,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecSF */
