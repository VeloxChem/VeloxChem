#ifndef ThreeCenterElectronRepulsionPrimRecFSP_hpp
#define ThreeCenterElectronRepulsionPrimRecFSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psp The primitive integrals buffer.
/// @param idx_eri_1_psp The primitive integrals buffer.
/// @param idx_eri_1_dss The primitive integrals buffer.
/// @param idx_eri_1_dsp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsp,
                                 size_t idx_eri_0_psp,
                                 size_t idx_eri_1_psp,
                                 size_t idx_eri_1_dss,
                                 size_t idx_eri_1_dsp,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSP_hpp */
