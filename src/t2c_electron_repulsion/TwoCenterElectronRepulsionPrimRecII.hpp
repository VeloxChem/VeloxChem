#ifndef TwoCenterElectronRepulsionPrimRecII
#define TwoCenterElectronRepulsionPrimRecII

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [I|1/|r-r'||I]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ii The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gi The index of integral in primitive integrals buffer.
/// @param idx_eri_1_gi The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_hi The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_ii(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ii,
                                const size_t idx_eri_0_gi,
                                const size_t idx_eri_1_gi,
                                const size_t idx_eri_1_hh,
                                const size_t idx_eri_1_hi,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecII */
