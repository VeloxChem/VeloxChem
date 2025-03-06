#ifndef ThreeCenterElectronRepulsionPrimRecHSD_hpp
#define ThreeCenterElectronRepulsionPrimRecHSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsd The primitive integrals buffer.
/// @param idx_eri_1_fsd The primitive integrals buffer.
/// @param idx_eri_1_gsp The primitive integrals buffer.
/// @param idx_eri_1_gsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsd,
                                 size_t idx_eri_0_fsd,
                                 size_t idx_eri_1_fsd,
                                 size_t idx_eri_1_gsp,
                                 size_t idx_eri_1_gsd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSD_hpp */
