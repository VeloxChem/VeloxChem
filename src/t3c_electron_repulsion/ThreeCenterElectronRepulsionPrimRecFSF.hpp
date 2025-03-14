#ifndef ThreeCenterElectronRepulsionPrimRecFSF_hpp
#define ThreeCenterElectronRepulsionPrimRecFSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psf The primitive integrals buffer.
/// @param idx_eri_1_psf The primitive integrals buffer.
/// @param idx_eri_1_dsd The primitive integrals buffer.
/// @param idx_eri_1_dsf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsf,
                                 size_t idx_eri_0_psf,
                                 size_t idx_eri_1_psf,
                                 size_t idx_eri_1_dsd,
                                 size_t idx_eri_1_dsf,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSF_hpp */
