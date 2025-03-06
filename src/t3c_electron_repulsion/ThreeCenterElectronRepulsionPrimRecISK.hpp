#ifndef ThreeCenterElectronRepulsionPrimRecISK_hpp
#define ThreeCenterElectronRepulsionPrimRecISK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_isk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsk The primitive integrals buffer.
/// @param idx_eri_1_gsk The primitive integrals buffer.
/// @param idx_eri_1_hsi The primitive integrals buffer.
/// @param idx_eri_1_hsk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_isk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isk,
                                 size_t idx_eri_0_gsk,
                                 size_t idx_eri_1_gsk,
                                 size_t idx_eri_1_hsi,
                                 size_t idx_eri_1_hsk,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISK_hpp */
