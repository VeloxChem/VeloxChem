#ifndef ThreeCenterElectronRepulsionPrimRecKSI_hpp
#define ThreeCenterElectronRepulsionPrimRecKSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [K|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ksi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_hsi The primitive integrals buffer.
/// @param idx_eri_1_hsi The primitive integrals buffer.
/// @param idx_eri_1_ish The primitive integrals buffer.
/// @param idx_eri_1_isi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ksi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksi,
                                 size_t idx_eri_0_hsi,
                                 size_t idx_eri_1_hsi,
                                 size_t idx_eri_1_ish,
                                 size_t idx_eri_1_isi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecKSI_hpp */
