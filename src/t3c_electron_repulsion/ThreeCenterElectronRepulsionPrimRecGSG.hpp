#ifndef ThreeCenterElectronRepulsionPrimRecGSG_hpp
#define ThreeCenterElectronRepulsionPrimRecGSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [G|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_gsg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_dsg The primitive integrals buffer.
/// @param idx_eri_1_dsg The primitive integrals buffer.
/// @param idx_eri_1_fsf The primitive integrals buffer.
/// @param idx_eri_1_fsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_gsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsg,
                                 size_t idx_eri_0_dsg,
                                 size_t idx_eri_1_dsg,
                                 size_t idx_eri_1_fsf,
                                 size_t idx_eri_1_fsg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecGSG_hpp */
