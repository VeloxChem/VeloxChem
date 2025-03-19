#ifndef ThreeCenterElectronRepulsionPrimRecGSI_hpp
#define ThreeCenterElectronRepulsionPrimRecGSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsi The primitive integrals buffer.
/// @param idx_eri_1_dsi The primitive integrals buffer.
/// @param idx_eri_1_fsh The primitive integrals buffer.
/// @param idx_eri_1_fsi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsi,
                                 size_t idx_eri_0_dsi,
                                 size_t idx_eri_1_dsi,
                                 size_t idx_eri_1_fsh,
                                 size_t idx_eri_1_fsi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSI_hpp */
