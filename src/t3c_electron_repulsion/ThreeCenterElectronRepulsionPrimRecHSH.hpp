#ifndef ThreeCenterElectronRepulsionPrimRecHSH_hpp
#define ThreeCenterElectronRepulsionPrimRecHSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsh The primitive integrals buffer.
/// @param idx_eri_1_fsh The primitive integrals buffer.
/// @param idx_eri_1_gsg The primitive integrals buffer.
/// @param idx_eri_1_gsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsh,
                                 size_t idx_eri_0_fsh,
                                 size_t idx_eri_1_fsh,
                                 size_t idx_eri_1_gsg,
                                 size_t idx_eri_1_gsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSH_hpp */
