#ifndef ElectronRepulsionGeom1010ContrRecGPXX_hpp
#define ElectronRepulsionGeom1010ContrRecGPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_gpxx The contracted integrals buffer.
/// @param idx_geom_00_fpxx The contracted integrals buffer.
/// @param idx_geom_10_fpxx The contracted integrals buffer.
/// @param idx_geom_10_fdxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_gpxx,
                                              const size_t idx_geom_0010_fpxx,
                                              const size_t idx_geom_1010_fpxx,
                                              const size_t idx_geom_1010_fdxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecGPXX_hpp */
