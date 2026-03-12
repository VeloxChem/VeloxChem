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
