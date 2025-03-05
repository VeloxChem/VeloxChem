#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXDD_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXDD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||DD)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xdd The contracted integrals buffer.
/// @param idx_xpd The contracted integrals buffer.
/// @param idx_geom_10_xpd The contracted integrals buffer.
/// @param idx_geom_10_xpf The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xdd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdd,
                                        const size_t idx_xpd,
                                        const size_t idx_geom_10_xpd,
                                        const size_t idx_geom_10_xpf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXDD_hpp */
