#ifndef ThreeCenterElectronRepulsionGeom010ContrRecXPH_hpp
#define ThreeCenterElectronRepulsionGeom010ContrRecXPH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dC^(1)(X|1/|r-r'||PH)  integral derivatives.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_xph The contracted integrals buffer.
/// @param idx_xsh The contracted integrals buffer.
/// @param idx_geom_10_xsh The contracted integrals buffer.
/// @param idx_geom_10_xsi The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_ket_geom010_electron_repulsion_xph(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xph,
                                        const size_t idx_xsh,
                                        const size_t idx_geom_10_xsh,
                                        const size_t idx_geom_10_xsi,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010ContrRecXPH_hpp */
