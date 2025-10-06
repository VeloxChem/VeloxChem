#ifndef ElectronRepulsionGeom0100ContrRecFKXX_hpp
#define ElectronRepulsionGeom0100ContrRecFKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_fkxx The contracted integrals buffer.
/// @param idx_dkxx The contracted integrals buffer.
/// @param idx_geom_01_dkxx The contracted integrals buffer.
/// @param idx_geom_01_dlxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_fkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_fkxx,
                                            const size_t idx_dkxx,
                                            const size_t idx_geom_01_dkxx,
                                            const size_t idx_geom_01_dlxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecFKXX_hpp */
