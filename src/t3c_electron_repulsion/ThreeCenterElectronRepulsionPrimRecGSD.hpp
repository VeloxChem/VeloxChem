#ifndef ThreeCenterElectronRepulsionPrimRecGSD_hpp
#define ThreeCenterElectronRepulsionPrimRecGSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsd The primitive integrals buffer.
/// @param idx_eri_1_dsd The primitive integrals buffer.
/// @param idx_eri_1_fsp The primitive integrals buffer.
/// @param idx_eri_1_fsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsd,
                                 size_t idx_eri_0_dsd,
                                 size_t idx_eri_1_dsd,
                                 size_t idx_eri_1_fsp,
                                 size_t idx_eri_1_fsd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSD_hpp */
