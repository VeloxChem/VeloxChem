#ifndef ElectronRepulsionGeom0100ContrRecHKXX_hpp
#define ElectronRepulsionGeom0100ContrRecHKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_hkxx The contracted integrals buffer.
/// @param idx_gkxx The contracted integrals buffer.
/// @param idx_geom_01_gkxx The contracted integrals buffer.
/// @param idx_geom_01_glxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_hkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hkxx,
                                            const size_t idx_gkxx,
                                            const size_t idx_geom_01_gkxx,
                                            const size_t idx_geom_01_glxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecHKXX_hpp */
