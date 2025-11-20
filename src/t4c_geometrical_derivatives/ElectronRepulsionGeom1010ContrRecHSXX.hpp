#ifndef ElectronRepulsionGeom1010ContrRecHSXX_hpp
#define ElectronRepulsionGeom1010ContrRecHSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HS|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_hsxx The contracted integrals buffer.
/// @param idx_geom_00_gsxx The contracted integrals buffer.
/// @param idx_geom_10_gsxx The contracted integrals buffer.
/// @param idx_geom_10_gpxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_hsxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_hsxx,
                                              const size_t idx_geom_0010_gsxx,
                                              const size_t idx_geom_1010_gsxx,
                                              const size_t idx_geom_1010_gpxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecHSXX_hpp */
