#ifndef ThreeCenterElectronRepulsionContrRecXDS_hpp
#define ThreeCenterElectronRepulsionContrRecXDS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||DS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xds The contracted integrals buffer.
/// @param idx_xps The contracted integrals buffer.
/// @param idx_xpp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xds(CSimdArray<double>& cbuffer,
                                const size_t idx_xds,
                                const size_t idx_xps,
                                const size_t idx_xpp,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXDS_hpp */
