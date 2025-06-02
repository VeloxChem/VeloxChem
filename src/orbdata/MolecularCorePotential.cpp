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

#include "MolecularCorePotential.hpp"

CMolecularCorePotential::CMolecularCorePotential()

    : _core_potentials{}

    , _indices{}

    , _atom_indices{}
{
}

CMolecularCorePotential::CMolecularCorePotential(const std::vector<CAtomCorePotential> &core_potentials, const std::vector<int> &indices, const std::vector<int> &atom_indices)
    
    : _core_potentials(core_potentials)

    , _indices(indices)

    , _atom_indices(atom_indices)
{
}

CMolecularCorePotential::CMolecularCorePotential(const CMolecularCorePotential &other)

    : _core_potentials(other._core_potentials)

    , _indices(other._indices)

    , _atom_indices(other._atom_indices)
{
}

CMolecularCorePotential::CMolecularCorePotential(CMolecularCorePotential &&other) noexcept

    : _core_potentials(std::move(other._core_potentials))

    , _indices(std::move(other._indices))

    , _atom_indices(std::move(other._atom_indices))
{
}

auto
CMolecularCorePotential::operator=(const CMolecularCorePotential &other) -> CMolecularCorePotential &
{
    _core_potentials = other._core_potentials;

    _indices = other._indices;
    
    _atom_indices = other._atom_indices; 

    return *this;
}

auto
CMolecularCorePotential::operator=(CMolecularCorePotential &&other) noexcept -> CMolecularCorePotential &
{
    if (this != &other)
    {
        _core_potentials = std::move(other._core_potentials);

        _indices = std::move(other._indices);
        
        _atom_indices = std::move(other._atom_indices);
    }

    return *this;
}

auto
CMolecularCorePotential::operator==(const CMolecularCorePotential &other) const -> bool
{
    if (_indices != other._indices)
    {
        return false;
    }
    else if (_atom_indices != other._atom_indices)
    {
        return false;
    }
    else
    {
        return _core_potentials == other._core_potentials;
    }
}

auto
CMolecularCorePotential::operator!=(const CMolecularCorePotential &other) const -> bool
{
    return !(*this == other);
}

auto
CMolecularCorePotential::add(const CAtomCorePotential &core_potential, const int iatom) -> void
{
    auto pos = std::ranges::find(_core_potentials, core_potential);

    if (pos == _core_potentials.end())
    {
        _indices.push_back(static_cast<int>(_core_potentials.size()));

        _core_potentials.push_back(core_potential);
    }
    else
    {
        _indices.push_back(static_cast<int>(std::distance(_core_potentials.begin(), pos)));
    }
    
    _atom_indices.push_back(iatom);
}

auto
CMolecularCorePotential::core_potentials() const -> std::vector<CAtomCorePotential>
{
    return _core_potentials;
}

auto
CMolecularCorePotential::indices() const -> std::vector<int>
{
    return _indices;
}

auto
CMolecularCorePotential::atomic_indices() const -> std::vector<int>
{
    return _atom_indices;
}

auto
CMolecularCorePotential::core_potential(const int index) const -> CAtomCorePotential
{
    return _core_potentials.at(index);
}
