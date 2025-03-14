#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXGD_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXGD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||GD)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xgd The contracted integrals buffer.
/// @param idx_xfd The contracted integrals buffer.
/// @param idx_geom_10_xfd The contracted integrals buffer.
/// @param idx_geom_10_xff The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xgd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xgd,
                                        const size_t idx_xfd,
                                        const size_t idx_geom_10_xfd,
                                        const size_t idx_geom_10_xff,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXGD_hpp */
