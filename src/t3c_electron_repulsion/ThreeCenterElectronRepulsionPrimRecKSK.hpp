#ifndef ThreeCenterElectronRepulsionPrimRecKSK_hpp
#define ThreeCenterElectronRepulsionPrimRecKSK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsk The primitive integrals buffer.
/// @param idx_eri_1_hsk The primitive integrals buffer.
/// @param idx_eri_1_isi The primitive integrals buffer.
/// @param idx_eri_1_isk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksk,
                                 size_t idx_eri_0_hsk,
                                 size_t idx_eri_1_hsk,
                                 size_t idx_eri_1_isi,
                                 size_t idx_eri_1_isk,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSK_hpp */
