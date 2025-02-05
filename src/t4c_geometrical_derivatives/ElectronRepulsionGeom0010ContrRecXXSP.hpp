#ifndef ElectronRepulsionGeom0010ContrRecXXSP_hpp
#define ElectronRepulsionGeom0010ContrRecXXSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||SP)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xxsp The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxsp The contracted integrals buffer.
/// @param idx_xxsd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_geom10_hrr_electron_repulsion_xxsp(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsp,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsp,
                                            const size_t idx_xxsd,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0010ContrRecXXSP_hpp */
