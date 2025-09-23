#ifndef ElectronRepulsionGeom0100ContrRecKGXX_hpp
#define ElectronRepulsionGeom0100ContrRecKGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_kgxx The contracted integrals buffer.
/// @param idx_igxx The contracted integrals buffer.
/// @param idx_geom_01_igxx The contracted integrals buffer.
/// @param idx_geom_01_ihxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_kgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kgxx,
                                            const size_t idx_igxx,
                                            const size_t idx_geom_01_igxx,
                                            const size_t idx_geom_01_ihxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKGXX_hpp */
