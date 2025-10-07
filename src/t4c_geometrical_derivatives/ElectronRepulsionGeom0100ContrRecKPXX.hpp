#ifndef ElectronRepulsionGeom0100ContrRecKPXX_hpp
#define ElectronRepulsionGeom0100ContrRecKPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_kpxx The contracted integrals buffer.
/// @param idx_ipxx The contracted integrals buffer.
/// @param idx_geom_01_ipxx The contracted integrals buffer.
/// @param idx_geom_01_idxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_kpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kpxx,
                                            const size_t idx_ipxx,
                                            const size_t idx_geom_01_ipxx,
                                            const size_t idx_geom_01_idxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKPXX_hpp */
