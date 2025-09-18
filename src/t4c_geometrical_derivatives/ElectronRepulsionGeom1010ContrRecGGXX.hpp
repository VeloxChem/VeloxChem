#ifndef ElectronRepulsionGeom1010ContrRecGGXX_hpp
#define ElectronRepulsionGeom1010ContrRecGGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_ggxx The contracted integrals buffer.
/// @param idx_geom_00_fgxx The contracted integrals buffer.
/// @param idx_geom_10_fgxx The contracted integrals buffer.
/// @param idx_geom_10_fhxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_ggxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_ggxx,
                                              const size_t idx_geom_0010_fgxx,
                                              const size_t idx_geom_1010_fgxx,
                                              const size_t idx_geom_1010_fhxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecGGXX_hpp */
