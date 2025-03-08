#ifndef ThreeCenterElectronRepulsionPrimRecKSD_hpp
#define ThreeCenterElectronRepulsionPrimRecKSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsd The primitive integrals buffer.
/// @param idx_eri_1_hsd The primitive integrals buffer.
/// @param idx_eri_1_isp The primitive integrals buffer.
/// @param idx_eri_1_isd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksd,
                                 size_t idx_eri_0_hsd,
                                 size_t idx_eri_1_hsd,
                                 size_t idx_eri_1_isp,
                                 size_t idx_eri_1_isd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSD_hpp */
