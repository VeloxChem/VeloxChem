#ifndef ThreeCenterElectronRepulsionPrimRecHSL_hpp
#define ThreeCenterElectronRepulsionPrimRecHSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [H|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_hsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_fsl The primitive integrals buffer.
/// @param idx_eri_1_fsl The primitive integrals buffer.
/// @param idx_eri_1_gsk The primitive integrals buffer.
/// @param idx_eri_1_gsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_hsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsl,
                                 size_t idx_eri_0_fsl,
                                 size_t idx_eri_1_fsl,
                                 size_t idx_eri_1_gsk,
                                 size_t idx_eri_1_gsl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecHSL_hpp */
