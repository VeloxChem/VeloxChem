#ifndef ThreeCenterElectronRepulsionContrRecXFS_hpp
#define ThreeCenterElectronRepulsionContrRecXFS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||FS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xfs The contracted integrals buffer.
/// @param idx_xds The contracted integrals buffer.
/// @param idx_xdp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xfs(CSimdArray<double>& cbuffer,
                                const size_t idx_xfs,
                                const size_t idx_xds,
                                const size_t idx_xdp,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXFS_hpp */
