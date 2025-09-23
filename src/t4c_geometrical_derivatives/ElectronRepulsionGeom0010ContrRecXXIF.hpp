#ifndef ElectronRepulsionGeom0010ContrRecXXIF_hpp
#define ElectronRepulsionGeom0010ContrRecXXIF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||IF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxif The contracted integrals buffer.
/// @param idx_xxhf The contracted integrals buffer.
/// @param idx_geom_10_xxhf The contracted integrals buffer.
/// @param idx_geom_10_xxhg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxif(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxif,
                                            const size_t idx_xxhf,
                                            const size_t idx_geom_10_xxhf,
                                            const size_t idx_geom_10_xxhg,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXIF_hpp */
