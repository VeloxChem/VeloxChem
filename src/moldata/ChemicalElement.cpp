//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "ChemicalElement.hpp"

#include <algorithm>
#include <iterator>

namespace chem_elem {

auto
name(const int id) -> std::string
{
    return _names.at(id);
}

auto
label(const int id) -> std::string
{
    auto label = _names.at(id);

    if (label.size() == 2)
    {
        label[1] = std::tolower(label[1]);
    }
    return label;
}

auto
identifier(const std::string &name) -> int
{
    if (auto it = std::ranges::find(_names, name); it != _names.end())
    {
        return static_cast<int>(std::distance(_names.begin(), it));
    }
    else
    {
        return -1;
    }
}

auto
mass(const int id) -> double
{
    return _masses.at(id);
}

auto
max_angular_momentum(const int id) -> int
{
    if ((id > 0) && (id < 5)) return 0;

    if ((id > 4) && (id < 21)) return 1;

    if ((id > 20) && (id < 57)) return 2;

    if ((id > 56) && (id < 87)) return 3;

    return -1;
}

auto
max_identifier() -> int
{
    return static_cast<int>(_names.size()) - 1;
}

}  // namespace chem_elem
