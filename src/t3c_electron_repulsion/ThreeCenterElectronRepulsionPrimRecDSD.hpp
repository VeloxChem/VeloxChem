#ifndef ThreeCenterElectronRepulsionPrimRecDSD_hpp
#define ThreeCenterElectronRepulsionPrimRecDSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssd The primitive integrals buffer.
/// @param idx_eri_1_ssd The primitive integrals buffer.
/// @param idx_eri_1_psp The primitive integrals buffer.
/// @param idx_eri_1_psd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsd,
                                 size_t idx_eri_0_ssd,
                                 size_t idx_eri_1_ssd,
                                 size_t idx_eri_1_psp,
                                 size_t idx_eri_1_psd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSD_hpp */
