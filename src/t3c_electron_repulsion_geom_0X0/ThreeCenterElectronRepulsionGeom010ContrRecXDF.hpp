#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXDF_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXDF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||DF)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xdf The contracted integrals buffer.
/// @param idx_xpf The contracted integrals buffer.
/// @param idx_geom_10_xpf The contracted integrals buffer.
/// @param idx_geom_10_xpg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xdf(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdf,
                                        const size_t idx_xpf,
                                        const size_t idx_geom_10_xpf,
                                        const size_t idx_geom_10_xpg,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXDF_hpp */
