#ifndef ThreeCenterElectronRepulsionPrimRecGSP_hpp
#define ThreeCenterElectronRepulsionPrimRecGSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsp The primitive integrals buffer.
/// @param idx_eri_1_dsp The primitive integrals buffer.
/// @param idx_eri_1_fss The primitive integrals buffer.
/// @param idx_eri_1_fsp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsp,
                                 size_t idx_eri_0_dsp,
                                 size_t idx_eri_1_dsp,
                                 size_t idx_eri_1_fss,
                                 size_t idx_eri_1_fsp,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSP_hpp */
