#ifndef ElectronRepulsionGeom0100ContrRecGPXX_hpp
#define ElectronRepulsionGeom0100ContrRecGPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_gpxx The contracted integrals buffer.
/// @param idx_fpxx The contracted integrals buffer.
/// @param idx_geom_01_fpxx The contracted integrals buffer.
/// @param idx_geom_01_fdxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gpxx,
                                            const size_t idx_fpxx,
                                            const size_t idx_geom_01_fpxx,
                                            const size_t idx_geom_01_fdxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecGPXX_hpp */
