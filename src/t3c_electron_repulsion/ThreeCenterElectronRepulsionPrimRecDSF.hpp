#ifndef ThreeCenterElectronRepulsionPrimRecDSF_hpp
#define ThreeCenterElectronRepulsionPrimRecDSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssf The primitive integrals buffer.
/// @param idx_eri_1_ssf The primitive integrals buffer.
/// @param idx_eri_1_psd The primitive integrals buffer.
/// @param idx_eri_1_psf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsf,
                                 size_t idx_eri_0_ssf,
                                 size_t idx_eri_1_ssf,
                                 size_t idx_eri_1_psd,
                                 size_t idx_eri_1_psf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSF_hpp */
