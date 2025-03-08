#ifndef ThreeCenterElectronRepulsionPrimRecGSF_hpp
#define ThreeCenterElectronRepulsionPrimRecGSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsf The primitive integrals buffer.
/// @param idx_eri_1_dsf The primitive integrals buffer.
/// @param idx_eri_1_fsd The primitive integrals buffer.
/// @param idx_eri_1_fsf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsf,
                                 size_t idx_eri_0_dsf,
                                 size_t idx_eri_1_dsf,
                                 size_t idx_eri_1_fsd,
                                 size_t idx_eri_1_fsf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSF_hpp */
