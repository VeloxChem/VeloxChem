#ifndef ElectronRepulsionGeom0100ContrRecHSXX_hpp
#define ElectronRepulsionGeom0100ContrRecHSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (HS|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_hsxx The contracted integrals buffer.
/// @param idx_gsxx The contracted integrals buffer.
/// @param idx_geom_01_gsxx The contracted integrals buffer.
/// @param idx_geom_01_gpxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_hsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hsxx,
                                            const size_t idx_gsxx,
                                            const size_t idx_geom_01_gsxx,
                                            const size_t idx_geom_01_gpxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecHSXX_hpp */
