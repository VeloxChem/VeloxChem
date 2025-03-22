#ifndef ThreeCenterElectronRepulsionPrimRecDSP_hpp
#define ThreeCenterElectronRepulsionPrimRecDSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssp The primitive integrals buffer.
/// @param idx_eri_1_ssp The primitive integrals buffer.
/// @param idx_eri_1_pss The primitive integrals buffer.
/// @param idx_eri_1_psp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsp,
                                 size_t idx_eri_0_ssp,
                                 size_t idx_eri_1_ssp,
                                 size_t idx_eri_1_pss,
                                 size_t idx_eri_1_psp,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSP_hpp */
