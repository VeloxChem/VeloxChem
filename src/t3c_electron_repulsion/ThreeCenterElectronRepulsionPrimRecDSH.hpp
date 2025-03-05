#ifndef ThreeCenterElectronRepulsionPrimRecDSH_hpp
#define ThreeCenterElectronRepulsionPrimRecDSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssh The primitive integrals buffer.
/// @param idx_eri_1_ssh The primitive integrals buffer.
/// @param idx_eri_1_psg The primitive integrals buffer.
/// @param idx_eri_1_psh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsh,
                                 size_t idx_eri_0_ssh,
                                 size_t idx_eri_1_ssh,
                                 size_t idx_eri_1_psg,
                                 size_t idx_eri_1_psh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSH_hpp */
