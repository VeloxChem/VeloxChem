#ifndef ThreeCenterElectronRepulsionGeom100ContrRecPXX_hpp
#define ThreeCenterElectronRepulsionGeom100ContrRecPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (P|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_100_pxx The contracted integrals buffer.
/// @param idx_sxx The contracted integrals buffer.
/// @param idx_dxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1_electron_repulsion_pxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_pxx,
                                      const size_t idx_sxx,
                                      const size_t idx_dxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100ContrRecPXX_hpp */
