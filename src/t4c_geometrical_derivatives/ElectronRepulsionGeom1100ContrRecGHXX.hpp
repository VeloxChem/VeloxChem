#ifndef ElectronRepulsionGeom1100ContrRecGHXX_hpp
#define ElectronRepulsionGeom1100ContrRecGHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GH|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_ghxx The contracted integrals buffer.
/// @param idx_geom_01_fhxx The contracted integrals buffer.
/// @param idx_geom_10_fhxx The contracted integrals buffer.
/// @param idx_geom_11_fhxx The contracted integrals buffer.
/// @param idx_geom_11_fixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_ghxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_ghxx,
                                            const size_t idx_geom_01_fhxx,
                                            const size_t idx_geom_10_fhxx,
                                            const size_t idx_geom_11_fhxx,
                                            const size_t idx_geom_11_fixx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecGHXX_hpp */
