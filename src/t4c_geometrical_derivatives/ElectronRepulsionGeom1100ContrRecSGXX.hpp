#ifndef ElectronRepulsionGeom1100ContrRecSGXX_hpp
#define ElectronRepulsionGeom1100ContrRecSGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_sgxx The contracted integrals buffer.
/// @param idx_sgxx The contracted integrals buffer.
/// @param idx_geom_01_sgxx The contracted integrals buffer.
/// @param idx_geom_01_shxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_sgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_sgxx,
                                            const size_t idx_sgxx,
                                            const size_t idx_geom_01_sgxx,
                                            const size_t idx_geom_01_shxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecSGXX_hpp */
