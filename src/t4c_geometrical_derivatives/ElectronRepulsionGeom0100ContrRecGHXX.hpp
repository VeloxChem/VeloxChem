#ifndef ElectronRepulsionGeom0100ContrRecGHXX_hpp
#define ElectronRepulsionGeom0100ContrRecGHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GH|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_ghxx The contracted integrals buffer.
/// @param idx_fhxx The contracted integrals buffer.
/// @param idx_geom_01_fhxx The contracted integrals buffer.
/// @param idx_geom_01_fixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_ghxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ghxx,
                                            const size_t idx_fhxx,
                                            const size_t idx_geom_01_fhxx,
                                            const size_t idx_geom_01_fixx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecGHXX_hpp */
