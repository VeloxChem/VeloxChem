#ifndef ThreeCenterElectronRepulsionPrimRecFSL_hpp
#define ThreeCenterElectronRepulsionPrimRecFSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psl The primitive integrals buffer.
/// @param idx_eri_1_psl The primitive integrals buffer.
/// @param idx_eri_1_dsk The primitive integrals buffer.
/// @param idx_eri_1_dsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsl,
                                 size_t idx_eri_0_psl,
                                 size_t idx_eri_1_psl,
                                 size_t idx_eri_1_dsk,
                                 size_t idx_eri_1_dsl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSL_hpp */
