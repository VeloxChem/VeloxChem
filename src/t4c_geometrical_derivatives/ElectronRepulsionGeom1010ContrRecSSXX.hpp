#ifndef ElectronRepulsionGeom1010ContrRecSSXX_hpp
#define ElectronRepulsionGeom1010ContrRecSSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SS|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_ssxx The contracted integrals buffer.
/// @param idx_geom_0010_ssxx The contracted integrals buffer.
/// @param idx_geom_0010_spxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_ssxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_ssxx,
                                              const size_t idx_geom_0010_ssxx,
                                              const size_t idx_geom_0010_spxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecSSXX_hpp */
