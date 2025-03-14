#ifndef ThreeCenterElectronRepulsionPrimRecHSI_hpp
#define ThreeCenterElectronRepulsionPrimRecHSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsi The primitive integrals buffer.
/// @param idx_eri_1_fsi The primitive integrals buffer.
/// @param idx_eri_1_gsh The primitive integrals buffer.
/// @param idx_eri_1_gsi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsi,
                                 size_t idx_eri_0_fsi,
                                 size_t idx_eri_1_fsi,
                                 size_t idx_eri_1_gsh,
                                 size_t idx_eri_1_gsi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSI_hpp */
