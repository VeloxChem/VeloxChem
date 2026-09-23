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


#include "EmbeddingRegion.hpp"

#include <algorithm>
#include <string>

#include "EmbeddingError.hpp"

CEmbeddingRegion::CEmbeddingRegion(const bool allows_polarizabilities)

    : _allows_polarizabilities(allows_polarizabilities)
{
}

auto
CEmbeddingRegion::allows_polarizabilities() const -> bool
{
    return _allows_polarizabilities;
}

auto
CEmbeddingRegion::add_force_field(const CPolarizableForceField &force_field) -> int
{
    embedding::require(_allows_polarizabilities || (!force_field.is_polarizable()),
                              std::string("EmbeddingRegion: The force field ") + force_field.get_name() +
                                  std::string(" carries a polarizability and this region takes none"));

    const auto index = index_of_force_field(force_field.get_name());

    if (index < 0)
    {
        _force_fields.push_back(force_field);

        _version++;

        return static_cast<int>(_force_fields.size()) - 1;
    }

    embedding::require(_force_fields[index].matches(force_field),
                              std::string("EmbeddingRegion: The region already holds a force field named ") +
                                  force_field.get_name() + std::string(" which carries other parameters"));

    return index;
}

auto
CEmbeddingRegion::add_molecule(const CMolecule &molecule, const int identifier) -> void
{
    embedding::require((identifier >= 0) && (identifier < number_of_force_fields()),
                              std::string("EmbeddingRegion: There is no force field of this identifier"));

    const auto nsites = static_cast<int>(_force_fields[identifier].number_of_sites());

    embedding::require(nsites == molecule.number_of_atoms(),
                              std::string("EmbeddingRegion: The force field ") + _force_fields[identifier].get_name() +
                                  std::string(" has one site for each of ") + std::to_string(nsites) +
                                  std::string(" atoms and this molecule has ") + std::to_string(molecule.number_of_atoms()));

    _molecules.push_back(molecule);

    _identifiers.push_back(identifier);

    _version++;
}

auto
CEmbeddingRegion::add_molecule(const CMolecule &molecule, const CPolarizableForceField &force_field) -> void
{
    add_molecule(molecule, add_force_field(force_field));
}

auto
CEmbeddingRegion::index_of_force_field(const std::string &name) const -> int
{
    for (size_t i = 0; i < _force_fields.size(); i++)
    {
        if (_force_fields[i].get_name() == name) return static_cast<int>(i);
    }

    return -1;
}

auto
CEmbeddingRegion::number_of_molecules() const -> int
{
    return static_cast<int>(_molecules.size());
}

auto
CEmbeddingRegion::number_of_force_fields() const -> int
{
    return static_cast<int>(_force_fields.size());
}

auto
CEmbeddingRegion::get_molecule(const int index) const -> const CMolecule &
{
    embedding::require((index >= 0) && (index < number_of_molecules()),
                              std::string("EmbeddingRegion: There is no molecule of this index"));

    return _molecules[index];
}

auto
CEmbeddingRegion::get_identifier(const int index) const -> int
{
    embedding::require((index >= 0) && (index < number_of_molecules()),
                              std::string("EmbeddingRegion: There is no molecule of this index"));

    return _identifiers[index];
}

auto
CEmbeddingRegion::get_force_field(const int index) const -> const CPolarizableForceField &
{
    embedding::require((index >= 0) && (index < number_of_force_fields()),
                              std::string("EmbeddingRegion: There is no force field of this index"));

    return _force_fields[index];
}

auto
CEmbeddingRegion::force_field_of(const int index) const -> const CPolarizableForceField &
{
    return _force_fields[get_identifier(index)];
}

auto
CEmbeddingRegion::get_molecules() const -> const std::vector<CMolecule> &
{
    return _molecules;
}

auto
CEmbeddingRegion::get_identifiers() const -> const std::vector<int> &
{
    return _identifiers;
}

auto
CEmbeddingRegion::get_force_fields() const -> const std::vector<CPolarizableForceField> &
{
    return _force_fields;
}

auto
CEmbeddingRegion::number_of_sites() const -> int
{
    int nsites = 0;

    for (const auto identifier : _identifiers) nsites += static_cast<int>(_force_fields[identifier].number_of_sites());

    return nsites;
}

auto
CEmbeddingRegion::number_of_polarizable_sites() const -> int
{
    std::vector<int> npolarizable(_force_fields.size(), 0);

    for (size_t i = 0; i < _force_fields.size(); i++)
    {
        for (const auto &site : _force_fields[i].get_sites())
        {
            if (site.is_polarizable()) npolarizable[i]++;
        }
    }

    int nsites = 0;

    for (const auto identifier : _identifiers) nsites += npolarizable[identifier];

    return nsites;
}

auto
CEmbeddingRegion::is_polarizable() const -> bool
{
    return number_of_polarizable_sites() > 0;
}

auto
CEmbeddingRegion::_carries(const CPolarizableSite &site, const int order) -> bool
{
    if (order == 0) return site.get_charge() != 0.0;

    if (site.get_order() < order) return false;

    if (order == 1)
    {
        const auto &dipole = site.get_dipole();

        return std::ranges::any_of(dipole, [](const double value) { return value != 0.0; });
    }

    const auto &quadrupole = site.get_quadrupole();

    return std::ranges::any_of(quadrupole, [](const double value) { return value != 0.0; });
}

auto
CEmbeddingRegion::_refuse_a_gaussian(const CPolarizableForceField &force_field,
                                     const size_t                  index,
                                     const CPolarizableSite       &site,
                                     const int                     order) -> void
{
    if (site.get_form() == pesite::point) return;

    // NOTE: the width is asked for only once the site is known to be a Gaussian
    // one. A point site has none and refuses the question, so a message which
    // names it cannot be written before the check.

    embedding::require(false,
                       std::string("EmbeddingRegion: The site ") + std::to_string(index + 1) +
                           std::string(" of the force field ") + force_field.get_name() +
                           std::string(" is a Gaussian of width ") + std::to_string(site.get_multipole_width()) +
                           std::string(", and the multipoles of a Gaussian are not the point multipoles of order ") +
                           std::to_string(order));
}

auto
CEmbeddingRegion::permanent_multipoles(const int order) const -> const TPermanentMultipoles &
{
    const auto ncomponents = multipoles::components(order);

    auto &gathered = _multipoles[static_cast<size_t>(order)];

    if (_has_multipoles[static_cast<size_t>(order)] && (_multipoles_version[static_cast<size_t>(order)] == _version))
    {
        return gathered;
    }

    gathered.order = order;

    gathered.coordinates.clear();

    gathered.values.clear();

    gathered.coordinates.reserve(3 * static_cast<size_t>(number_of_sites()));

    gathered.values.reserve(ncomponents * static_cast<size_t>(number_of_sites()));

    for (size_t imol = 0; imol < _molecules.size(); imol++)
    {
        const auto &force_field = _force_fields[static_cast<size_t>(_identifiers[imol])];

        const auto &coordinates = _molecules[imol].coordinates();

        for (size_t isite = 0; isite < force_field.number_of_sites(); isite++)
        {
            const auto &site = force_field.get_site(isite);

            if (!_carries(site, order)) continue;

            _refuse_a_gaussian(force_field, isite, site, order);

            const auto xyz = coordinates[isite].coordinates();

            gathered.coordinates.insert(gathered.coordinates.end(), xyz.begin(), xyz.end());

            if (order == 0)
            {
                gathered.values.push_back(site.get_charge());
            }
            else if (order == 1)
            {
                const auto &dipole = site.get_dipole();

                gathered.values.insert(gathered.values.end(), dipole.begin(), dipole.end());
            }
            else
            {
                const auto &quadrupole = site.get_quadrupole();

                gathered.values.insert(gathered.values.end(), quadrupole.begin(), quadrupole.end());
            }
        }
    }

    _multipoles_version[static_cast<size_t>(order)] = _version;

    _has_multipoles[static_cast<size_t>(order)] = true;

    return gathered;
}

auto
CEmbeddingRegion::version() const -> size_t
{
    return _version;
}
