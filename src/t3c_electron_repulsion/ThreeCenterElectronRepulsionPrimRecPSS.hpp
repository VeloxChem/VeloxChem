#ifndef ThreeCenterElectronRepulsionPrimRecPSS_hpp
#define ThreeCenterElectronRepulsionPrimRecPSS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [P|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_pss The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
auto
comp_prim_electron_repulsion_pss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_pss,
                                 size_t idx_eri_1_sss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecPSS_hpp */
