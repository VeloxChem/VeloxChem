#ifndef ElectronRepulsionGeom0010ContrRecXXDP_hpp
#define ElectronRepulsionGeom0010ContrRecXXDP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||DP)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxdp The contracted integrals buffer.
/// @param idx_xxpp The contracted integrals buffer.
/// @param idx_geom_10_xxpp The contracted integrals buffer.
/// @param idx_geom_10_xxpd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxdp(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxdp,
                                            const size_t idx_xxpp,
                                            const size_t idx_geom_10_xxpp,
                                            const size_t idx_geom_10_xxpd,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXDP_hpp */
