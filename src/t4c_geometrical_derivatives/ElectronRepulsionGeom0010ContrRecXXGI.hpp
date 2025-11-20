#ifndef ElectronRepulsionGeom0010ContrRecXXGI_hpp
#define ElectronRepulsionGeom0010ContrRecXXGI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||GI)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxgi The contracted integrals buffer.
/// @param idx_xxfi The contracted integrals buffer.
/// @param idx_geom_10_xxfi The contracted integrals buffer.
/// @param idx_geom_10_xxfk The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxgi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgi,
                                            const size_t idx_xxfi,
                                            const size_t idx_geom_10_xxfi,
                                            const size_t idx_geom_10_xxfk,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXGI_hpp */
