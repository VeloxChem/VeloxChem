#ifndef ThreeCenterElectronRepulsionPrimRecHSF_hpp
#define ThreeCenterElectronRepulsionPrimRecHSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsf The primitive integrals buffer.
/// @param idx_eri_1_fsf The primitive integrals buffer.
/// @param idx_eri_1_gsd The primitive integrals buffer.
/// @param idx_eri_1_gsf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsf,
                                 size_t idx_eri_0_fsf,
                                 size_t idx_eri_1_fsf,
                                 size_t idx_eri_1_gsd,
                                 size_t idx_eri_1_gsf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSF_hpp */
