#ifndef AngularMomentum_hpp
#define AngularMomentum_hpp

#include <cstdint>
#include <string>

#include "T2Pair.hpp"

namespace angmom {  // angmom namespace

/**
 Determines number of spherical components for given angular momentum.

 @param angmom the angular momentum.
 @return the number of spherical components.
 */
inline auto
to_SphericalComponents(const int64_t angmom) -> int64_t
{
    return 2 * angmom + 1;
}

/**
 Determines number of spherical components for given pair of angular momentums.

  @param bra_angmom the angular momentum on bra side.
  @param ket_angmom the angular momentum on ket side.
  @return the number of spherical components.
*/
inline auto
to_SphericalComponents(const int64_t bra_angmom, const int64_t ket_angmom) -> int64_t
{
    return (2 * bra_angmom + 1) * (2 * ket_angmom + 1);
}

/**
 Determines number of Cartesian components for given angular momentum.

 @param angmom the angular momentum.
 @return the number of Cartesian momentum.
 */
inline auto
to_CartesianComponents(const int64_t angmom) -> int64_t
{
    return (angmom + 1) * (angmom + 2) / 2;
}

/**
 Determines number of Cartesian components for given pair of angular momentums.

 @param bra_angmom the angular momentum on bra side.
 @param ket_angmom the angular momentum on ket side.
 @return the number of Cartesian momentum.
 */
inline auto
to_CartesianComponents(const int64_t bra_angmom, const int64_t ket_angmom) -> int64_t
{
    return ((bra_angmom + 1) * (bra_angmom + 2) / 2) * ((ket_angmom + 1) * (ket_angmom + 2) / 2);
}

/**
 Gets string representation of spherical angular momentum component.

 @param angmom the spherical angular momentum.
 @param component the spherical component of angular momentum.
 @return the string of angular momentum component.
 */
auto getStringOfAngularMomentum(const int64_t angmom, const int64_t component) -> std::string;

}  // namespace angmom

#endif /* AngularMomentum_hpp */
