#ifndef ElectronRepulsionGeom2000ContrRecDIXX_hpp
#define ElectronRepulsionGeom2000ContrRecDIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DI|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_20_dixx The contracted integrals buffer.
/// @param idx_geom_10_pixx The contracted integrals buffer.
/// @param idx_geom_20_pixx The contracted integrals buffer.
/// @param idx_geom_20_pkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom20_hrr_electron_repulsion_dixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dixx,
                                            const size_t idx_geom_10_pixx,
                                            const size_t idx_geom_20_pixx,
                                            const size_t idx_geom_20_pkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom2000ContrRecDIXX_hpp */
