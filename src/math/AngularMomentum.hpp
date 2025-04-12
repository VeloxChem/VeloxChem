//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
