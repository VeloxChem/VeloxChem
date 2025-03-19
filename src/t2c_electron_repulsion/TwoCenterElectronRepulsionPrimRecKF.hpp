#ifndef TwoCenterElectronRepulsionPrimRecKF
#define TwoCenterElectronRepulsionPrimRecKF

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [K|1/|r-r'||F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_kf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hf The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hf The index of integral in primitive integrals buffer.
/// @param idx_eri_1_id The index of integral in primitive integrals buffer.
/// @param idx_eri_1_if The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_kf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kf,
                                const size_t idx_eri_0_hf,
                                const size_t idx_eri_1_hf,
                                const size_t idx_eri_1_id,
                                const size_t idx_eri_1_if,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecKF */
