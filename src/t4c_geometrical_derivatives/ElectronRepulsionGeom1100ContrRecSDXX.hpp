#ifndef ElectronRepulsionGeom1100ContrRecSDXX_hpp
#define ElectronRepulsionGeom1100ContrRecSDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SD|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_sdxx The contracted integrals buffer.
/// @param idx_sdxx The contracted integrals buffer.
/// @param idx_geom_01_sdxx The contracted integrals buffer.
/// @param idx_geom_01_sfxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_sdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_sdxx,
                                            const size_t idx_sdxx,
                                            const size_t idx_geom_01_sdxx,
                                            const size_t idx_geom_01_sfxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecSDXX_hpp */
