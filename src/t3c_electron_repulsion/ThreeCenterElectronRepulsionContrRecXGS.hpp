#ifndef ThreeCenterElectronRepulsionContrRecXGS_hpp
#define ThreeCenterElectronRepulsionContrRecXGS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||GS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xgs The contracted integrals buffer.
/// @param idx_xfs The contracted integrals buffer.
/// @param idx_xfp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xgs(CSimdArray<double>& cbuffer,
                                const size_t idx_xgs,
                                const size_t idx_xfs,
                                const size_t idx_xfp,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXGS_hpp */
