#ifndef ThreeCenterElectronRepulsionPrimRecISD_hpp
#define ThreeCenterElectronRepulsionPrimRecISD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [I|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_isd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_gsd The primitive integrals buffer.
/// @param idx_eri_1_gsd The primitive integrals buffer.
/// @param idx_eri_1_hsp The primitive integrals buffer.
/// @param idx_eri_1_hsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_isd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isd,
                                 size_t idx_eri_0_gsd,
                                 size_t idx_eri_1_gsd,
                                 size_t idx_eri_1_hsp,
                                 size_t idx_eri_1_hsd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecISD_hpp */
