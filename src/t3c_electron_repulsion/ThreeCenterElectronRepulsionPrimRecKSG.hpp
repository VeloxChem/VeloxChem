#ifndef ThreeCenterElectronRepulsionPrimRecKSG_hpp
#define ThreeCenterElectronRepulsionPrimRecKSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsg The primitive integrals buffer.
/// @param idx_eri_1_hsg The primitive integrals buffer.
/// @param idx_eri_1_isf The primitive integrals buffer.
/// @param idx_eri_1_isg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksg,
                                 size_t idx_eri_0_hsg,
                                 size_t idx_eri_1_hsg,
                                 size_t idx_eri_1_isf,
                                 size_t idx_eri_1_isg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSG_hpp */
