#ifndef ElectronRepulsionGeom1010ContrRecHHXX_hpp
#define ElectronRepulsionGeom1010ContrRecHHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HH|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_hhxx The contracted integrals buffer.
/// @param idx_geom_00_ghxx The contracted integrals buffer.
/// @param idx_geom_10_ghxx The contracted integrals buffer.
/// @param idx_geom_10_gixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_hhxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_hhxx,
                                              const size_t idx_geom_0010_ghxx,
                                              const size_t idx_geom_1010_ghxx,
                                              const size_t idx_geom_1010_gixx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecHHXX_hpp */
