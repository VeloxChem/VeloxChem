#ifndef ElectronRepulsionGeom0100ContrRecKKXX_hpp
#define ElectronRepulsionGeom0100ContrRecKKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_kkxx The contracted integrals buffer.
/// @param idx_ikxx The contracted integrals buffer.
/// @param idx_geom_01_ikxx The contracted integrals buffer.
/// @param idx_geom_01_ilxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_kkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kkxx,
                                            const size_t idx_ikxx,
                                            const size_t idx_geom_01_ikxx,
                                            const size_t idx_geom_01_ilxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKKXX_hpp */
