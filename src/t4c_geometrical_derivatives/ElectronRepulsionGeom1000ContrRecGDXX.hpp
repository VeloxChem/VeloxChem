#ifndef ElectronRepulsionGeom1000ContrRecGDXX_hpp
#define ElectronRepulsionGeom1000ContrRecGDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GD|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_gdxx The contracted integrals buffer.
/// @param idx_fdxx The contracted integrals buffer.
/// @param idx_geom_10_fdxx The contracted integrals buffer.
/// @param idx_geom_10_ffxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_gdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gdxx,
                                            const size_t idx_fdxx,
                                            const size_t idx_geom_10_fdxx,
                                            const size_t idx_geom_10_ffxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecGDXX_hpp */
