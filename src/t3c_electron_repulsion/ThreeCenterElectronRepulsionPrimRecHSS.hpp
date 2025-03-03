#ifndef ThreeCenterElectronRepulsionPrimRecHSS_hpp
#define ThreeCenterElectronRepulsionPrimRecHSS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hss The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fss The primitive integrals buffer.
/// @param idx_eri_1_fss The primitive integrals buffer.
/// @param idx_eri_1_gss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hss,
                                 size_t idx_eri_0_fss,
                                 size_t idx_eri_1_fss,
                                 size_t idx_eri_1_gss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSS_hpp */
