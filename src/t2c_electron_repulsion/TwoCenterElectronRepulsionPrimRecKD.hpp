#ifndef TwoCenterElectronRepulsionPrimRecKD
#define TwoCenterElectronRepulsionPrimRecKD

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [K|1/|r-r'||D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_kd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hd The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hd The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ip The index of integral in primitive integrals buffer.
/// @param idx_eri_1_id The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_kd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kd,
                                const size_t idx_eri_0_hd,
                                const size_t idx_eri_1_hd,
                                const size_t idx_eri_1_ip,
                                const size_t idx_eri_1_id,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecKD */
