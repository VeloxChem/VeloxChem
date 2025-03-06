#ifndef ThreeCenterElectronRepulsionPrimRecISH_hpp
#define ThreeCenterElectronRepulsionPrimRecISH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ish The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsh The primitive integrals buffer.
/// @param idx_eri_1_gsh The primitive integrals buffer.
/// @param idx_eri_1_hsg The primitive integrals buffer.
/// @param idx_eri_1_hsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ish(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ish,
                                 size_t idx_eri_0_gsh,
                                 size_t idx_eri_1_gsh,
                                 size_t idx_eri_1_hsg,
                                 size_t idx_eri_1_hsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISH_hpp */
