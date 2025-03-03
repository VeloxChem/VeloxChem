#ifndef ThreeCenterElectronRepulsionPrimRecKSH_hpp
#define ThreeCenterElectronRepulsionPrimRecKSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsh The primitive integrals buffer.
/// @param idx_eri_1_hsh The primitive integrals buffer.
/// @param idx_eri_1_isg The primitive integrals buffer.
/// @param idx_eri_1_ish The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksh,
                                 size_t idx_eri_0_hsh,
                                 size_t idx_eri_1_hsh,
                                 size_t idx_eri_1_isg,
                                 size_t idx_eri_1_ish,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSH_hpp */
