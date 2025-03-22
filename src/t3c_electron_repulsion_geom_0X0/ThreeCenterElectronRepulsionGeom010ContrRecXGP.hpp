#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXGP_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXGP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||GP)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xgp The contracted integrals buffer.
/// @param idx_xfp The contracted integrals buffer.
/// @param idx_geom_10_xfp The contracted integrals buffer.
/// @param idx_geom_10_xfd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xgp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xgp,
                                        const size_t idx_xfp,
                                        const size_t idx_geom_10_xfp,
                                        const size_t idx_geom_10_xfd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXGP_hpp */
