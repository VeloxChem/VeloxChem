#ifndef ThreeCenterElectronRepulsionPrimRecPSF_hpp
#define ThreeCenterElectronRepulsionPrimRecPSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [P|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_psf The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ssd The primitive integrals buffer.
/// @param idx_eri_1_ssf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_psf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psf,
                                 size_t idx_eri_1_ssd,
                                 size_t idx_eri_1_ssf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecPSF_hpp */
