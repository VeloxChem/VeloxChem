#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXDH_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXDH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||DH)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xdh The contracted integrals buffer.
/// @param idx_xph The contracted integrals buffer.
/// @param idx_geom_10_xph The contracted integrals buffer.
/// @param idx_geom_10_xpi The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xdh(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdh,
                                        const size_t idx_xph,
                                        const size_t idx_geom_10_xph,
                                        const size_t idx_geom_10_xpi,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXDH_hpp */
