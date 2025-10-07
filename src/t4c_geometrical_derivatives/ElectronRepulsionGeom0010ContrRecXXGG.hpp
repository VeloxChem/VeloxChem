#ifndef ElectronRepulsionGeom0010ContrRecXXGG_hpp
#define ElectronRepulsionGeom0010ContrRecXXGG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||GG)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxgg The contracted integrals buffer.
/// @param idx_xxfg The contracted integrals buffer.
/// @param idx_geom_10_xxfg The contracted integrals buffer.
/// @param idx_geom_10_xxfh The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxgg(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgg,
                                            const size_t idx_xxfg,
                                            const size_t idx_geom_10_xxfg,
                                            const size_t idx_geom_10_xxfh,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXGG_hpp */
