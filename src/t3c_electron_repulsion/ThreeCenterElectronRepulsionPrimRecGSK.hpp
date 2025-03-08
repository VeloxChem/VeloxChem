#ifndef ThreeCenterElectronRepulsionPrimRecGSK_hpp
#define ThreeCenterElectronRepulsionPrimRecGSK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsk The primitive integrals buffer.
/// @param idx_eri_1_dsk The primitive integrals buffer.
/// @param idx_eri_1_fsi The primitive integrals buffer.
/// @param idx_eri_1_fsk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsk,
                                 size_t idx_eri_0_dsk,
                                 size_t idx_eri_1_dsk,
                                 size_t idx_eri_1_fsi,
                                 size_t idx_eri_1_fsk,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSK_hpp */
