#ifndef ElectronRepulsionGeom0100ContrRecPKXX_hpp
#define ElectronRepulsionGeom0100ContrRecPKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_pkxx The contracted integrals buffer.
/// @param idx_skxx The contracted integrals buffer.
/// @param idx_geom_01_skxx The contracted integrals buffer.
/// @param idx_geom_01_slxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_pkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_pkxx,
                                            const size_t idx_skxx,
                                            const size_t idx_geom_01_skxx,
                                            const size_t idx_geom_01_slxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecPKXX_hpp */
