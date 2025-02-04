#ifndef TwoCenterElectronRepulsionPrimRecGF
#define TwoCenterElectronRepulsionPrimRecGF

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [G|1/|r-r'||F]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_df The index of integral in primitive integrals buffer.
/// @param idx_eri_1_df The index of integral in primitive integrals buffer.
/// @param idx_eri_1_fd The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ff The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_gf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gf,
                                const size_t idx_eri_0_df,
                                const size_t idx_eri_1_df,
                                const size_t idx_eri_1_fd,
                                const size_t idx_eri_1_ff,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecGF */
