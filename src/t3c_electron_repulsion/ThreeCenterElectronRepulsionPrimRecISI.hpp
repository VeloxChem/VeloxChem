#ifndef ThreeCenterElectronRepulsionPrimRecISI_hpp
#define ThreeCenterElectronRepulsionPrimRecISI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_isi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsi The primitive integrals buffer.
/// @param idx_eri_1_gsi The primitive integrals buffer.
/// @param idx_eri_1_hsh The primitive integrals buffer.
/// @param idx_eri_1_hsi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_isi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isi,
                                 size_t idx_eri_0_gsi,
                                 size_t idx_eri_1_gsi,
                                 size_t idx_eri_1_hsh,
                                 size_t idx_eri_1_hsi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISI_hpp */
