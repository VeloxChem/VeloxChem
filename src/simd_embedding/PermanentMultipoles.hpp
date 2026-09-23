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


#ifndef PermanentMultipoles_hpp
#define PermanentMultipoles_hpp

#include <cstddef>
#include <vector>

namespace multipoles {  // multipoles namespace

/// @brief Gets how many numbers a multipole of an order is written with.
/// @param order The order: zero for a charge, one for a dipole, two for a
/// quadrupole.
/// @return One, three or six.
auto components(const int order) -> size_t;

}  // namespace multipoles

/// @brief Struct TPermanentMultipoles holds the permanent multipoles of one
/// order of an environment, put onto the molecules which carry them.
///
/// @note This is what a force field and an instance of it make together. The
/// force fields carry no coordinates, because one of them serves every molecule
/// of its kind; the molecules carry no parameters. Walking the two together,
/// site against atom and in that order, is what puts a multipole somewhere, and
/// this is the result of that walk.
///
/// @note Flat arrays, because that is what the integral drivers take: three
/// coordinates for each site and one, three or six values.
///
/// @note A site which carries nothing of this order is not here. An absent
/// moment and a zero one are different things to whoever contracts them, and
/// leaving the absent ones out is also most of the array for an environment of
/// water, whose hydrogens carry a charge and nothing else.
struct TPermanentMultipoles
{
    /// @brief The order of these multipoles.
    int order = 0;

    /// @brief The position of each site, three coordinates each, in bohr.
    std::vector<double> coordinates;

    /// @brief The multipole of each site, one, three or six numbers each.
    std::vector<double> values;

    /// @brief Gets the number of sites.
    auto number_of_sites() const -> size_t;

    /// @brief Gets how many numbers each site's multipole is written with.
    auto components() const -> size_t;
};

#endif /* PermanentMultipoles_hpp */
