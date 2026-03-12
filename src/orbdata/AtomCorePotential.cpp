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

#include "AtomCorePotential.hpp"

CAtomCorePotential::CAtomCorePotential()
    
    : _local_potential{}

    , _projected_potentials{}

    , _angular_momentums{}

    , _core_electrons{0}
{
}

CAtomCorePotential::CAtomCorePotential(const CBaseCorePotential              &local_potential,
                                       const std::vector<CBaseCorePotential> &projected_potentials,
                                       const std::vector<int>                &angular_momentums,
                                       const int                             core_electrons)

    : _local_potential(local_potential)

    , _projected_potentials(projected_potentials)

    , _angular_momentums(angular_momentums)

    , _core_electrons(core_electrons)
{
}

CAtomCorePotential::CAtomCorePotential(const CAtomCorePotential &other)

    : _local_potential(other._local_potential)

    , _projected_potentials(other._projected_potentials)

    , _angular_momentums(other._angular_momentums)

    , _core_electrons(other._core_electrons)
{
}

CAtomCorePotential::CAtomCorePotential(CAtomCorePotential &&other) noexcept

    : _local_potential{}

    , _projected_potentials{}

    , _angular_momentums{}

    , _core_electrons{0}
{
    std::swap(_local_potential, other._local_potential);

    std::swap(_projected_potentials, other._projected_potentials);

    std::swap(_angular_momentums, other._angular_momentums);
    
    std::swap(_core_electrons, other._core_electrons);
}

auto
CAtomCorePotential::operator=(const CAtomCorePotential &other) -> CAtomCorePotential &
{
    _local_potential = other._local_potential;

    _projected_potentials = other._projected_potentials;

    _angular_momentums = other._angular_momentums;
    
    _core_electrons = other._core_electrons;

    return *this;
}

auto
CAtomCorePotential::operator=(CAtomCorePotential &&other) noexcept -> CAtomCorePotential &
{
    std::swap(_local_potential, other._local_potential);

    std::swap(_projected_potentials, other._projected_potentials);

    std::swap(_angular_momentums, other._angular_momentums);

    std::swap(_core_electrons, other._core_electrons);
    
    return *this;
}

auto
CAtomCorePotential::operator==(const CAtomCorePotential &other) const -> bool
{
    if (_core_electrons != other._core_electrons)
    {
        return false;
    }
    else if (_angular_momentums != other._angular_momentums)
    {
        return false;
    }
    else if (_local_potential != other._local_potential)
    {
        return false;
    }
    else
    {
        return _projected_potentials == other._projected_potentials;
    }
}

auto
CAtomCorePotential::operator!=(const CAtomCorePotential &other) const -> bool
{
    return !(*this == other);
}

auto
CAtomCorePotential::set_local_potential(const CBaseCorePotential &local_potential) -> void
{
    _local_potential = local_potential;
}

auto
CAtomCorePotential::set_projected_potentials(const std::vector<CBaseCorePotential> &projected_potentials,
                                             const std::vector<int>                &angular_momentums) -> void
{
    _projected_potentials = projected_potentials;
    
    _angular_momentums = angular_momentums;
}

auto
CAtomCorePotential::add_projected_potential(const CBaseCorePotential &projected_potential,
                                            const int                 angular_momentum) -> void
{
    _projected_potentials.push_back(projected_potential);
    
    _angular_momentums.push_back(angular_momentum);
}

auto
CAtomCorePotential::set_number_core_electrons(const int core_electrons) -> void
{
    _core_electrons = core_electrons;
}

auto
CAtomCorePotential::get_local_potential() const -> CBaseCorePotential
{
    return _local_potential;
}

auto
CAtomCorePotential::get_projected_potentials() const -> std::vector<CBaseCorePotential>
{
    return _projected_potentials;
}

auto
CAtomCorePotential::get_angular_momentums() const -> std::vector<int>
{
    return _angular_momentums;
}

auto
CAtomCorePotential::number_of_core_electrons() const -> int
{
    return _core_electrons;
}
