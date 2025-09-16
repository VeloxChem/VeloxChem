#ifndef ElectronRepulsionGeom0100ContrRecDIXX_hpp
#define ElectronRepulsionGeom0100ContrRecDIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DI|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_dixx The contracted integrals buffer.
/// @param idx_pixx The contracted integrals buffer.
/// @param idx_geom_01_pixx The contracted integrals buffer.
/// @param idx_geom_01_pkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_dixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dixx,
                                            const size_t idx_pixx,
                                            const size_t idx_geom_01_pixx,
                                            const size_t idx_geom_01_pkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecDIXX_hpp */
