#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXFD_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXFD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||FD)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xfd The contracted integrals buffer.
/// @param idx_xdd The contracted integrals buffer.
/// @param idx_geom_10_xdd The contracted integrals buffer.
/// @param idx_geom_10_xdf The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xfd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfd,
                                        const size_t idx_xdd,
                                        const size_t idx_geom_10_xdd,
                                        const size_t idx_geom_10_xdf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXFD_hpp */
