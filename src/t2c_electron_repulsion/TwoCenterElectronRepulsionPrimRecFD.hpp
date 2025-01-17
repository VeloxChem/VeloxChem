#ifndef TwoCenterElectronRepulsionPrimRecFD
#define TwoCenterElectronRepulsionPrimRecFD

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [F|1/|r-r'||D]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_pd The index of integral in primitive integrals buffer.
/// @param idx_eri_1_pd The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dp The index of integral in primitive integrals buffer.
/// @param idx_eri_1_dd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_fd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fd,
                                const size_t idx_eri_0_pd,
                                const size_t idx_eri_1_pd,
                                const size_t idx_eri_1_dp,
                                const size_t idx_eri_1_dd,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecFD */