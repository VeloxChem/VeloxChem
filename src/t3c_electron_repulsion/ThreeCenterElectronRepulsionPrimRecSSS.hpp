#ifndef ThreeCenterElectronRepulsionPrimRecSSP_hpp
#define ThreeCenterElectronRepulsionPrimRecSSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [S|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_ovl The index of combined overlap factors.
/// @param bf_data The Boys function data.
/// @param idx_bvals The index of Boys function data.
auto comp_prim_electron_repulsion_sss(CSimdArray<double>&       pbuffer,
                                       const size_t              idx_eri_0_sss,
                                       CSimdArray<double>&       factors,
                                       const size_t              idx_ovl,
                                       const CSimdArray<double>& bf_data,
                                       const size_t              idx_bvals) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecSSP_hpp */
