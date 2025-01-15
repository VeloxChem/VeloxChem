#ifndef ElectronRepulsionGeom0100ContrRecDPXX_hpp
#define ElectronRepulsionGeom0100ContrRecDPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_dpxx The contracted integrals buffer.
/// @param idx_ppxx The contracted integrals buffer.
/// @param idx_geom_01_ppxx The contracted integrals buffer.
/// @param idx_geom_01_pdxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_dpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dpxx,
                                            const size_t idx_ppxx,
                                            const size_t idx_geom_01_ppxx,
                                            const size_t idx_geom_01_pdxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecDPXX_hpp */
