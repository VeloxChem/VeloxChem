#ifndef ThreeCenterElectronRepulsionPrimRecFSD_hpp
#define ThreeCenterElectronRepulsionPrimRecFSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psd The primitive integrals buffer.
/// @param idx_eri_1_psd The primitive integrals buffer.
/// @param idx_eri_1_dsp The primitive integrals buffer.
/// @param idx_eri_1_dsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsd,
                                 size_t idx_eri_0_psd,
                                 size_t idx_eri_1_psd,
                                 size_t idx_eri_1_dsp,
                                 size_t idx_eri_1_dsd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSD_hpp */
