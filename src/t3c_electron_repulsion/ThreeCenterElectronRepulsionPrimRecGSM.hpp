#ifndef ThreeCenterElectronRepulsionPrimRecGSM_hpp
#define ThreeCenterElectronRepulsionPrimRecGSM_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SM]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsm The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsm The primitive integrals buffer.
/// @param idx_eri_1_dsm The primitive integrals buffer.
/// @param idx_eri_1_fsl The primitive integrals buffer.
/// @param idx_eri_1_fsm The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsm,
                                 size_t idx_eri_0_dsm,
                                 size_t idx_eri_1_dsm,
                                 size_t idx_eri_1_fsl,
                                 size_t idx_eri_1_fsm,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSM_hpp */
