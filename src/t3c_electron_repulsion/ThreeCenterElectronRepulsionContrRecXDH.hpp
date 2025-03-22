#ifndef ThreeCenterElectronRepulsionContrRecXDH_hpp
#define ThreeCenterElectronRepulsionContrRecXDH_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||DH)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xdh The contracted integrals buffer.
/// @param idx_xph The contracted integrals buffer.
/// @param idx_xpi The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xdh(CSimdArray<double>& cbuffer,
                                const size_t idx_xdh,
                                const size_t idx_xph,
                                const size_t idx_xpi,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXDH_hpp */
