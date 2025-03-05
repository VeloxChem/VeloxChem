#ifndef ThreeCenterElectronRepulsionPrimRecSSF_hpp
#define ThreeCenterElectronRepulsionPrimRecSSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [S|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ssf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssp The primitive integrals buffer.
/// @param idx_eri_1_ssp The primitive integrals buffer.
/// @param idx_eri_0_ssd The primitive integrals buffer.
/// @param idx_eri_1_ssd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ssf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssf,
                                 size_t idx_eri_0_ssp,
                                 size_t idx_eri_1_ssp,
                                 size_t idx_eri_0_ssd,
                                 size_t idx_eri_1_ssd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecSSF_hpp */
