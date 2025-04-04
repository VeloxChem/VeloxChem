#ifndef ThreeCenterElectronRepulsionPrimRecFSM_hpp
#define ThreeCenterElectronRepulsionPrimRecFSM_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SM]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsm The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psm The primitive integrals buffer.
/// @param idx_eri_1_psm The primitive integrals buffer.
/// @param idx_eri_1_dsl The primitive integrals buffer.
/// @param idx_eri_1_dsm The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsm,
                                 size_t idx_eri_0_psm,
                                 size_t idx_eri_1_psm,
                                 size_t idx_eri_1_dsl,
                                 size_t idx_eri_1_dsm,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSM_hpp */
