#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXFS_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXFS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||FS)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xfs The contracted integrals buffer.
/// @param idx_xds The contracted integrals buffer.
/// @param idx_geom_10_xds The contracted integrals buffer.
/// @param idx_geom_10_xdp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xfs(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfs,
                                        const size_t idx_xds,
                                        const size_t idx_geom_10_xds,
                                        const size_t idx_geom_10_xdp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXFS_hpp */
