#ifndef ElectronRepulsionGeom2000ContrRecHHXX_hpp
#define ElectronRepulsionGeom2000ContrRecHHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HH|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_20_hhxx The contracted integrals buffer.
/// @param idx_geom_10_ghxx The contracted integrals buffer.
/// @param idx_geom_20_ghxx The contracted integrals buffer.
/// @param idx_geom_20_gixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom20_hrr_electron_repulsion_hhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_hhxx,
                                            const size_t idx_geom_10_ghxx,
                                            const size_t idx_geom_20_ghxx,
                                            const size_t idx_geom_20_gixx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom2000ContrRecHHXX_hpp */
