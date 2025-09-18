#ifndef ElectronRepulsionGeom1010ContrRecHGXX_hpp
#define ElectronRepulsionGeom1010ContrRecHGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_hgxx The contracted integrals buffer.
/// @param idx_geom_00_ggxx The contracted integrals buffer.
/// @param idx_geom_10_ggxx The contracted integrals buffer.
/// @param idx_geom_10_ghxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_hgxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_hgxx,
                                              const size_t idx_geom_0010_ggxx,
                                              const size_t idx_geom_1010_ggxx,
                                              const size_t idx_geom_1010_ghxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecHGXX_hpp */
