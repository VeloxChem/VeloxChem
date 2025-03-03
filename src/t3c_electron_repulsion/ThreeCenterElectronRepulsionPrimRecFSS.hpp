#ifndef ThreeCenterElectronRepulsionPrimRecFSS_hpp
#define ThreeCenterElectronRepulsionPrimRecFSS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fss The index of integral in primitive integrals buffer.
/// @param idx_eri_0_pss The primitive integrals buffer.
/// @param idx_eri_1_pss The primitive integrals buffer.
/// @param idx_eri_1_dss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fss,
                                 size_t idx_eri_0_pss,
                                 size_t idx_eri_1_pss,
                                 size_t idx_eri_1_dss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSS_hpp */
