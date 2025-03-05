#ifndef ThreeCenterElectronRepulsionPrimRecGSL_hpp
#define ThreeCenterElectronRepulsionPrimRecGSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsl The primitive integrals buffer.
/// @param idx_eri_1_dsl The primitive integrals buffer.
/// @param idx_eri_1_fsk The primitive integrals buffer.
/// @param idx_eri_1_fsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsl,
                                 size_t idx_eri_0_dsl,
                                 size_t idx_eri_1_dsl,
                                 size_t idx_eri_1_fsk,
                                 size_t idx_eri_1_fsl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSL_hpp */
