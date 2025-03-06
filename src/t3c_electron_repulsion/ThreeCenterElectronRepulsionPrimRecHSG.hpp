#ifndef ThreeCenterElectronRepulsionPrimRecHSG_hpp
#define ThreeCenterElectronRepulsionPrimRecHSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsg The primitive integrals buffer.
/// @param idx_eri_1_fsg The primitive integrals buffer.
/// @param idx_eri_1_gsf The primitive integrals buffer.
/// @param idx_eri_1_gsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsg,
                                 size_t idx_eri_0_fsg,
                                 size_t idx_eri_1_fsg,
                                 size_t idx_eri_1_gsf,
                                 size_t idx_eri_1_gsg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSG_hpp */
