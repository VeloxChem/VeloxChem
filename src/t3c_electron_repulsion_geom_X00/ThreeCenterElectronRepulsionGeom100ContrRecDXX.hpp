#ifndef ThreeCenterElectronRepulsionGeom100ContrRecDXX_hpp
#define ThreeCenterElectronRepulsionGeom100ContrRecDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (D|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_100_dxx The contracted integrals buffer.
/// @param idx_pxx The contracted integrals buffer.
/// @param idx_fxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1_electron_repulsion_dxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_dxx,
                                      const size_t idx_pxx,
                                      const size_t idx_fxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100ContrRecDXX_hpp */
