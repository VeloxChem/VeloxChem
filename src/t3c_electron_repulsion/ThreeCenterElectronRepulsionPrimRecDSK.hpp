#ifndef ThreeCenterElectronRepulsionPrimRecDSK_hpp
#define ThreeCenterElectronRepulsionPrimRecDSK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssk The primitive integrals buffer.
/// @param idx_eri_1_ssk The primitive integrals buffer.
/// @param idx_eri_1_psi The primitive integrals buffer.
/// @param idx_eri_1_psk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsk,
                                 size_t idx_eri_0_ssk,
                                 size_t idx_eri_1_ssk,
                                 size_t idx_eri_1_psi,
                                 size_t idx_eri_1_psk,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSK_hpp */
