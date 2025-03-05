#ifndef ThreeCenterElectronRepulsionPrimRecKSF_hpp
#define ThreeCenterElectronRepulsionPrimRecKSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsf The primitive integrals buffer.
/// @param idx_eri_1_hsf The primitive integrals buffer.
/// @param idx_eri_1_isd The primitive integrals buffer.
/// @param idx_eri_1_isf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksf,
                                 size_t idx_eri_0_hsf,
                                 size_t idx_eri_1_hsf,
                                 size_t idx_eri_1_isd,
                                 size_t idx_eri_1_isf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSF_hpp */
