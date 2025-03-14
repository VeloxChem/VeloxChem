#ifndef ThreeCenterElectronRepulsionPrimRecHSK_hpp
#define ThreeCenterElectronRepulsionPrimRecHSK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsk The primitive integrals buffer.
/// @param idx_eri_1_fsk The primitive integrals buffer.
/// @param idx_eri_1_gsi The primitive integrals buffer.
/// @param idx_eri_1_gsk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsk,
                                 size_t idx_eri_0_fsk,
                                 size_t idx_eri_1_fsk,
                                 size_t idx_eri_1_gsi,
                                 size_t idx_eri_1_gsk,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSK_hpp */
