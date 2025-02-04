#ifndef ThreeCenterElectronRepulsionPrimRecDSL_hpp
#define ThreeCenterElectronRepulsionPrimRecDSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssl The primitive integrals buffer.
/// @param idx_eri_1_ssl The primitive integrals buffer.
/// @param idx_eri_1_psk The primitive integrals buffer.
/// @param idx_eri_1_psl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsl,
                                 size_t idx_eri_0_ssl,
                                 size_t idx_eri_1_ssl,
                                 size_t idx_eri_1_psk,
                                 size_t idx_eri_1_psl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSL_hpp */
