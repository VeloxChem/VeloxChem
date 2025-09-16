#ifndef ElectronRepulsionGeom1100ContrRecFDXX_hpp
#define ElectronRepulsionGeom1100ContrRecFDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FD|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_fdxx The contracted integrals buffer.
/// @param idx_geom_01_ddxx The contracted integrals buffer.
/// @param idx_geom_10_ddxx The contracted integrals buffer.
/// @param idx_geom_11_ddxx The contracted integrals buffer.
/// @param idx_geom_11_dfxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_fdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_fdxx,
                                            const size_t idx_geom_01_ddxx,
                                            const size_t idx_geom_10_ddxx,
                                            const size_t idx_geom_11_ddxx,
                                            const size_t idx_geom_11_dfxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecFDXX_hpp */
