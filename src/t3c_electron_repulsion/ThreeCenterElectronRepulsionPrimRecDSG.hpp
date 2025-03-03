#ifndef ThreeCenterElectronRepulsionPrimRecDSG_hpp
#define ThreeCenterElectronRepulsionPrimRecDSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [D|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_dsg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssg The primitive integrals buffer.
/// @param idx_eri_1_ssg The primitive integrals buffer.
/// @param idx_eri_1_psf The primitive integrals buffer.
/// @param idx_eri_1_psg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_dsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsg,
                                 size_t idx_eri_0_ssg,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_1_psf,
                                 size_t idx_eri_1_psg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecDSG_hpp */
