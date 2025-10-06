#ifndef ElectronRepulsionGeom0100ContrRecFGXX_hpp
#define ElectronRepulsionGeom0100ContrRecFGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_fgxx The contracted integrals buffer.
/// @param idx_dgxx The contracted integrals buffer.
/// @param idx_geom_01_dgxx The contracted integrals buffer.
/// @param idx_geom_01_dhxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_fgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_fgxx,
                                            const size_t idx_dgxx,
                                            const size_t idx_geom_01_dgxx,
                                            const size_t idx_geom_01_dhxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecFGXX_hpp */
