#ifndef ThreeCenterElectronRepulsionPrimRecSSP_hpp
#define ThreeCenterElectronRepulsionPrimRecSSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [S|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ssp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sss The primitive integrals buffer.
/// @param idx_eri_1_sss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
auto
comp_prim_electron_repulsion_ssp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssp,
                                 size_t idx_eri_0_sss,
                                 size_t idx_eri_1_sss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecSSP_hpp */
