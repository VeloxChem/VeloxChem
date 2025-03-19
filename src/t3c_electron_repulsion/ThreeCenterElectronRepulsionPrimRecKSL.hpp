#ifndef ThreeCenterElectronRepulsionPrimRecKSL_hpp
#define ThreeCenterElectronRepulsionPrimRecKSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsl The primitive integrals buffer.
/// @param idx_eri_1_hsl The primitive integrals buffer.
/// @param idx_eri_1_isk The primitive integrals buffer.
/// @param idx_eri_1_isl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksl,
                                 size_t idx_eri_0_hsl,
                                 size_t idx_eri_1_hsl,
                                 size_t idx_eri_1_isk,
                                 size_t idx_eri_1_isl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSL_hpp */
