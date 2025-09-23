#ifndef ElectronRepulsionGeom0010ContrRecXXDF_hpp
#define ElectronRepulsionGeom0010ContrRecXXDF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||DF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxdf The contracted integrals buffer.
/// @param idx_xxpf The contracted integrals buffer.
/// @param idx_geom_10_xxpf The contracted integrals buffer.
/// @param idx_geom_10_xxpg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxdf(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxdf,
                                            const size_t idx_xxpf,
                                            const size_t idx_geom_10_xxpf,
                                            const size_t idx_geom_10_xxpg,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXDF_hpp */
