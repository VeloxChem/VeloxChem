#ifndef ElectronRepulsionGeom0100ContrRecKDXX_hpp
#define ElectronRepulsionGeom0100ContrRecKDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KD|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_kdxx The contracted integrals buffer.
/// @param idx_idxx The contracted integrals buffer.
/// @param idx_geom_01_idxx The contracted integrals buffer.
/// @param idx_geom_01_ifxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_kdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kdxx,
                                            const size_t idx_idxx,
                                            const size_t idx_geom_01_idxx,
                                            const size_t idx_geom_01_ifxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKDXX_hpp */
