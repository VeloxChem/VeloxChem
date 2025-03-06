#ifndef ThreeCenterElectronRepulsionPrimRecPSH_hpp
#define ThreeCenterElectronRepulsionPrimRecPSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [P|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_psh The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ssg The primitive integrals buffer.
/// @param idx_eri_1_ssh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_psh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psh,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_1_ssh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecPSH_hpp */
