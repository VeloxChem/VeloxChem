#ifndef ElectronRepulsionGeom1000ContrRecGIXX_hpp
#define ElectronRepulsionGeom1000ContrRecGIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GI|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_gixx The contracted integrals buffer.
/// @param idx_fixx The contracted integrals buffer.
/// @param idx_geom_10_fixx The contracted integrals buffer.
/// @param idx_geom_10_fkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_gixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gixx,
                                            const size_t idx_fixx,
                                            const size_t idx_geom_10_fixx,
                                            const size_t idx_geom_10_fkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecGIXX_hpp */
