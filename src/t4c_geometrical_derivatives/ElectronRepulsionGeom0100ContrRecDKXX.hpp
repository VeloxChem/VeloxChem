#ifndef ElectronRepulsionGeom0100ContrRecDKXX_hpp
#define ElectronRepulsionGeom0100ContrRecDKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_dkxx The contracted integrals buffer.
/// @param idx_pkxx The contracted integrals buffer.
/// @param idx_geom_01_pkxx The contracted integrals buffer.
/// @param idx_geom_01_plxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_dkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dkxx,
                                            const size_t idx_pkxx,
                                            const size_t idx_geom_01_pkxx,
                                            const size_t idx_geom_01_plxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecDKXX_hpp */
