#ifndef ThreeCenterElectronRepulsionPrimRecISG_hpp
#define ThreeCenterElectronRepulsionPrimRecISG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_isg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsg The primitive integrals buffer.
/// @param idx_eri_1_gsg The primitive integrals buffer.
/// @param idx_eri_1_hsf The primitive integrals buffer.
/// @param idx_eri_1_hsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_isg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isg,
                                 size_t idx_eri_0_gsg,
                                 size_t idx_eri_1_gsg,
                                 size_t idx_eri_1_hsf,
                                 size_t idx_eri_1_hsg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISG_hpp */
