#ifndef ElectronRepulsionGeom1100ContrRecPGXX_hpp
#define ElectronRepulsionGeom1100ContrRecPGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_pgxx The contracted integrals buffer.
/// @param idx_geom_01_sgxx The contracted integrals buffer.
/// @param idx_geom_10_sgxx The contracted integrals buffer.
/// @param idx_geom_11_sgxx The contracted integrals buffer.
/// @param idx_geom_11_shxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_pgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_pgxx,
                                            const size_t idx_geom_01_sgxx,
                                            const size_t idx_geom_10_sgxx,
                                            const size_t idx_geom_11_sgxx,
                                            const size_t idx_geom_11_shxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecPGXX_hpp */
