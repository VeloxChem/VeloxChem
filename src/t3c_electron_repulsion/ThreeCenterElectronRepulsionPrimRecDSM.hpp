#ifndef ThreeCenterElectronRepulsionPrimRecDSM_hpp
#define ThreeCenterElectronRepulsionPrimRecDSM_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SM]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsm The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssm The primitive integrals buffer.
/// @param idx_eri_1_ssm The primitive integrals buffer.
/// @param idx_eri_1_psl The primitive integrals buffer.
/// @param idx_eri_1_psm The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsm,
                                 size_t idx_eri_0_ssm,
                                 size_t idx_eri_1_ssm,
                                 size_t idx_eri_1_psl,
                                 size_t idx_eri_1_psm,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSM_hpp */
