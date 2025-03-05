#ifndef ThreeCenterElectronRepulsionGeom100ContrRecHXX_hpp
#define ThreeCenterElectronRepulsionGeom100ContrRecHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (H|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_100_hxx The contracted integrals buffer.
/// @param idx_gxx The contracted integrals buffer.
/// @param idx_ixx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1_electron_repulsion_hxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_hxx,
                                      const size_t idx_gxx,
                                      const size_t idx_ixx,
                                      const int c_angmom,
                                      const int d_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100ContrRecHXX_hpp */
