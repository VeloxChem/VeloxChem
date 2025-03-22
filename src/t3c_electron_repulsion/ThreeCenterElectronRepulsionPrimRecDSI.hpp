#ifndef ThreeCenterElectronRepulsionPrimRecDSI_hpp
#define ThreeCenterElectronRepulsionPrimRecDSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssi The primitive integrals buffer.
/// @param idx_eri_1_ssi The primitive integrals buffer.
/// @param idx_eri_1_psh The primitive integrals buffer.
/// @param idx_eri_1_psi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsi,
                                 size_t idx_eri_0_ssi,
                                 size_t idx_eri_1_ssi,
                                 size_t idx_eri_1_psh,
                                 size_t idx_eri_1_psi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSI_hpp */
