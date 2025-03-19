#ifndef ThreeCenterElectronRepulsionGeom100ContrRecIXX_hpp
#define ThreeCenterElectronRepulsionGeom100ContrRecIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (I|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_100_ixx The contracted integrals buffer.
/// @param idx_hxx The contracted integrals buffer.
/// @param idx_kxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1_electron_repulsion_ixx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_ixx,
                                      const size_t idx_hxx,
                                      const size_t idx_kxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100ContrRecIXX_hpp */
