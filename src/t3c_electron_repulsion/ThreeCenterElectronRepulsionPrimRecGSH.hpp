#ifndef ThreeCenterElectronRepulsionPrimRecGSH_hpp
#define ThreeCenterElectronRepulsionPrimRecGSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsh The primitive integrals buffer.
/// @param idx_eri_1_dsh The primitive integrals buffer.
/// @param idx_eri_1_fsg The primitive integrals buffer.
/// @param idx_eri_1_fsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsh,
                                 size_t idx_eri_0_dsh,
                                 size_t idx_eri_1_dsh,
                                 size_t idx_eri_1_fsg,
                                 size_t idx_eri_1_fsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSH_hpp */
