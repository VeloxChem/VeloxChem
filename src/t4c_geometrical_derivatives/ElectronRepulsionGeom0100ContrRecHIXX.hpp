#ifndef ElectronRepulsionGeom0100ContrRecHIXX_hpp
#define ElectronRepulsionGeom0100ContrRecHIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HI|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_hixx The contracted integrals buffer.
/// @param idx_gixx The contracted integrals buffer.
/// @param idx_geom_01_gixx The contracted integrals buffer.
/// @param idx_geom_01_gkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_hixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hixx,
                                            const size_t idx_gixx,
                                            const size_t idx_geom_01_gixx,
                                            const size_t idx_geom_01_gkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecHIXX_hpp */
