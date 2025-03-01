#ifndef ThreeCenterElectronRepulsionPrimRecISF_hpp
#define ThreeCenterElectronRepulsionPrimRecISF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_isf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsf The primitive integrals buffer.
/// @param idx_eri_1_gsf The primitive integrals buffer.
/// @param idx_eri_1_hsd The primitive integrals buffer.
/// @param idx_eri_1_hsf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_isf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isf,
                                 size_t idx_eri_0_gsf,
                                 size_t idx_eri_1_gsf,
                                 size_t idx_eri_1_hsd,
                                 size_t idx_eri_1_hsf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISF_hpp */
