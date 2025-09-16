#ifndef ElectronRepulsionGeom0100ContrRecKHXX_hpp
#define ElectronRepulsionGeom0100ContrRecKHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KH|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_khxx The contracted integrals buffer.
/// @param idx_ihxx The contracted integrals buffer.
/// @param idx_geom_01_ihxx The contracted integrals buffer.
/// @param idx_geom_01_iixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_khxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_khxx,
                                            const size_t idx_ihxx,
                                            const size_t idx_geom_01_ihxx,
                                            const size_t idx_geom_01_iixx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKHXX_hpp */
