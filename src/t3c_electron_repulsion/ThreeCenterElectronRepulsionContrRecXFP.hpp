#ifndef ThreeCenterElectronRepulsionContrRecXFP_hpp
#define ThreeCenterElectronRepulsionContrRecXFP_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||FP)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xfp The contracted integrals buffer.
/// @param idx_xdp The contracted integrals buffer.
/// @param idx_xdd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xfp(CSimdArray<double>& cbuffer,
                                const size_t idx_xfp,
                                const size_t idx_xdp,
                                const size_t idx_xdd,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXFP_hpp */
