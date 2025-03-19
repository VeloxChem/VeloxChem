#ifndef ThreeCenterElectronRepulsionGeom100ContrRecGXX_hpp
#define ThreeCenterElectronRepulsionGeom100ContrRecGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (G|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_100_gxx The contracted integrals buffer.
/// @param idx_fxx The contracted integrals buffer.
/// @param idx_hxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1_electron_repulsion_gxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_gxx,
                                      const size_t idx_fxx,
                                      const size_t idx_hxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100ContrRecGXX_hpp */
