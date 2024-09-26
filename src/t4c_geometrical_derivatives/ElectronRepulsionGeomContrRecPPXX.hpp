#ifndef ElectronRepulsionGeomContrRecPPXX_hpp
#define ElectronRepulsionGeomContrRecPPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (PP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_ppxx The contracted integrals buffer.
/// @param idx_geom_spxx The contracted integrals buffer.
/// @param idx_geom_sdxx The contracted integrals buffer.
/// @param idx_spxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto comp_bra_geom_hrr_electron_repulsion_ppxx(CSimdArray<double>&   cbuffer,
                                               const size_t          idx_geom_ppxx,
                                               const size_t          idx_geom_spxx,
                                               const size_t          idx_geom_sdxx,
                                               const size_t          idx_spxx,
                                               const TPoint<double>& r_ab,
                                               const int             c_angmom,
                                               const int             d_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionGeomContrRecPPXX_hpp */