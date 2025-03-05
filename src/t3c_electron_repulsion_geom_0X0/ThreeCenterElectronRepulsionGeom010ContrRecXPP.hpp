#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXPP_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXPP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||PP)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xpp The contracted integrals buffer.
/// @param idx_xsp The contracted integrals buffer.
/// @param idx_geom_10_xsp The contracted integrals buffer.
/// @param idx_geom_10_xsd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xpp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xpp,
                                        const size_t idx_xsp,
                                        const size_t idx_geom_10_xsp,
                                        const size_t idx_geom_10_xsd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXPP_hpp */
