#ifndef ThreeCenterElectronRepulsionContrRecXFG_hpp
#define ThreeCenterElectronRepulsionContrRecXFG_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||FG)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xfg The contracted integrals buffer.
/// @param idx_xdg The contracted integrals buffer.
/// @param idx_xdh The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xfg(CSimdArray<double>& cbuffer,
                                const size_t idx_xfg,
                                const size_t idx_xdg,
                                const size_t idx_xdh,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXFG_hpp */
