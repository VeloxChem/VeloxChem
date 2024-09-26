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

#include "Molecule.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <ranges>

#include "AtomicRadii.hpp"
#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "CustomViews.hpp"
#include "StringFormat.hpp"

CMolecule::CMolecule()

    : _charge{0.0}

    , _multiplicity{1}

    , _coordinates{}

    , _identifiers{}

    , _atom_basis_labels{}
{
}

CMolecule::CMolecule(const std::vector<int> &identifiers, const std::vector<TPoint<double>> &coordinates, const std::string &unit)

    : _charge{0.0}

    , _multiplicity{1}

    , _coordinates(coordinates)

    , _identifiers(identifiers)
{
    for (std::vector<int>::size_type i = 0; i < _identifiers.size(); i++)
    {
        _atom_basis_labels.push_back(std::string(""));
    }

    if (_is_angstrom(unit))
    {
        std::ranges::for_each(_coordinates, [=](TPoint<double> &pnt) { pnt.scale(1.0 / units::bohr_in_angstrom()); });
    }
}

CMolecule::CMolecule(const std::vector<int> &identifiers, const std::vector<TPoint<double>> &coordinates, const std::string &unit, const std::vector<std::string>& atom_basis_labels)

    : _charge{0.0}

    , _multiplicity{1}

    , _coordinates(coordinates)

    , _identifiers(identifiers)

    , _atom_basis_labels(atom_basis_labels)
{
    if (_is_angstrom(unit))
    {
        std::ranges::for_each(_coordinates, [=](TPoint<double> &pnt) { pnt.scale(1.0 / units::bohr_in_angstrom()); });
    }
}

CMolecule::CMolecule(const CMolecule &molecule_one, const CMolecule &molecule_two)
{
    _identifiers.reserve(molecule_one._identifiers.size() + molecule_two._identifiers.size());

    std::ranges::copy(molecule_one._identifiers, std::back_inserter(_identifiers));

    std::ranges::copy(molecule_two._identifiers, std::back_inserter(_identifiers));

    _coordinates.reserve(molecule_one._coordinates.size() + molecule_two._coordinates.size());

    std::ranges::copy(molecule_one._coordinates, std::back_inserter(_coordinates));

    std::ranges::copy(molecule_two._coordinates, std::back_inserter(_coordinates));

    _atom_basis_labels.reserve(molecule_one._atom_basis_labels.size() + molecule_two._atom_basis_labels.size());

    std::ranges::copy(molecule_one._atom_basis_labels, std::back_inserter(_atom_basis_labels));

    std::ranges::copy(molecule_two._atom_basis_labels, std::back_inserter(_atom_basis_labels));

    _charge = molecule_one._charge + molecule_two._charge;

    const auto spin_a = (molecule_one._multiplicity - 1) / 2;

    const auto spin_b = (molecule_two._multiplicity - 1) / 2;

    _multiplicity = 2 * (spin_a + spin_b) + 1;
}

CMolecule::CMolecule(const CMolecule &other)

    : _charge(other._charge)

    , _multiplicity(other._multiplicity)

    , _coordinates(other._coordinates)

    , _identifiers(other._identifiers)

    , _atom_basis_labels(other._atom_basis_labels)
{
}

CMolecule::CMolecule(CMolecule &&other) noexcept

    : _charge{0.0}

    , _multiplicity{1}

    , _coordinates{}

    , _identifiers{}

    , _atom_basis_labels{}
{
    std::swap(_charge, other._charge);

    std::swap(_multiplicity, other._multiplicity);

    std::swap(_coordinates, other._coordinates);

    std::swap(_identifiers, other._identifiers);

    std::swap(_atom_basis_labels, other._atom_basis_labels);
}

auto
CMolecule::operator=(const CMolecule &other) -> CMolecule &
{
    _charge = other._charge;

    _multiplicity = other._multiplicity;

    _coordinates = other._coordinates;

    _identifiers = other._identifiers;

    _atom_basis_labels = other._atom_basis_labels;

    return *this;
}

auto
CMolecule::operator=(CMolecule &&other) noexcept -> CMolecule &
{
    std::swap(_charge, other._charge);

    std::swap(_multiplicity, other._multiplicity);

    std::swap(_coordinates, other._coordinates);

    std::swap(_identifiers, other._identifiers);

    std::swap(_atom_basis_labels, other._atom_basis_labels);

    return *this;
}

auto
CMolecule::operator==(const CMolecule &other) const -> bool
{
    if (_multiplicity != other._multiplicity)
    {
        return false;
    }
    else if (!mathfunc::equal(_charge, other._charge, 1.0e-12, 1.0e-12))
    {
        return false;
    }
    else if (_identifiers != other._identifiers)
    {
        return false;
    }
    else if (_atom_basis_labels != other._atom_basis_labels)
    {
        return false;
    }
    else
    {
        return _coordinates == other._coordinates;
    }
}

auto
CMolecule::add_atom(const int identifier, const TPoint<double>& coordinates, const std::string& unit) -> void
{
    add_atom(identifier, coordinates, unit, std::string(""));
}

auto
CMolecule::add_atom(const int identifier, const TPoint<double> &coordinates, const std::string &unit, const std::string& atom_basis_label) -> void
{
    _identifiers.push_back(identifier);

    _atom_basis_labels.push_back(atom_basis_label);

    _coordinates.push_back(coordinates);

    if (_is_angstrom(unit))
    {
        _coordinates.back().scale(1.0 / units::bohr_in_angstrom());
    }
}

auto
CMolecule::slice(const std::vector<int> &atoms) const -> CMolecule
{
    CMolecule molfrag;

    std::ranges::for_each(atoms, [&](const int i) { molfrag.add_atom(_identifiers.at(i), _coordinates.at(i), "au", _atom_basis_labels.at(i)); });

    if ((molfrag.number_of_electrons() % 2) == 1) molfrag.set_multiplicity(2);

    return molfrag;
}

auto
CMolecule::set_charge(const double charge) -> void
{
    _charge = charge;
}

auto
CMolecule::set_multiplicity(const int multiplicity) -> void
{
    _multiplicity = multiplicity;
}

auto
CMolecule::get_charge() const -> double
{
    return _charge;
}

auto
CMolecule::get_multiplicity() const -> int
{
    return _multiplicity;
}

auto
CMolecule::number_of_atoms() const -> int
{
    return static_cast<int>(_identifiers.size());
}

auto
CMolecule::number_of_atoms(const int identifier) const -> int
{
    return static_cast<int>(std::ranges::count(_identifiers, identifier));
}

auto
CMolecule::number_of_atoms(const int iatom, const int natoms, const int identifier) const -> int
{
    int count = 0;

    std::ranges::for_each(std::views::iota(iatom, iatom + natoms), [&](const int i) {
        if (identifier == _identifiers.at(i)) count++;
    });

    return count;
}

auto
CMolecule::elemental_composition() const -> std::set<int>
{
    return std::set<int>(_identifiers.begin(), _identifiers.end());
}

auto
CMolecule::number_of_electrons() const -> int
{
    return std::accumulate(_identifiers.begin(), _identifiers.end(), -static_cast<int>(_charge));
}

auto
CMolecule::identifiers() const -> std::vector<int>
{
    return _identifiers;
}

auto
CMolecule::atom_basis_labels() const -> std::vector<std::string>
{
    return _atom_basis_labels;
}

auto
CMolecule::coordinates(const std::string &unit) const -> std::vector<TPoint<double>>
{
    if (_is_angstrom(unit))
    {
        std::vector<TPoint<double>> coords;

        coords.reserve(_coordinates.size());

        auto fact = units::bohr_in_angstrom();

        std::ranges::transform(_coordinates, std::back_inserter(coords), [&](const TPoint<double> &pnt) {
            TPoint<double> rpnt = pnt;
            rpnt.scale(fact);
            return rpnt;
        });

        return coords;
    }
    else
    {
        return _coordinates;
    }
}

auto
CMolecule::charges() const -> std::vector<double>
{
    std::vector<double> charges;

    charges.reserve(_identifiers.size());

    std::ranges::transform(_identifiers, std::back_inserter(charges), [](const int i) { return static_cast<double>(i); });

    return charges;
}

auto
CMolecule::masses() const -> std::vector<double>
{
    std::vector<double> masses;

    masses.reserve(_identifiers.size());

    std::ranges::transform(_identifiers, std::back_inserter(masses), [](const int i) { return chem_elem::mass(i); });

    return masses;
}

auto
CMolecule::labels() const -> std::vector<std::string>
{
    std::vector<std::string> labels;

    labels.reserve(_identifiers.size());

    std::ranges::transform(_identifiers, std::back_inserter(labels), [](const int i) { return chem_elem::label(i); });

    return labels;
}

auto
CMolecule::label(const int iatom) const -> std::string
{
    return chem_elem::label(_identifiers.at(iatom));
}

auto
CMolecule::atom_coordinates(const int iatom, const std::string &unit) const -> TPoint<double>
{
    if (_is_angstrom(unit))
    {
        auto rpnt = _coordinates.at(iatom);

        rpnt.scale(units::bohr_in_angstrom());

        return rpnt;
    }
    else
    {
        return _coordinates.at(iatom);
    }
}

auto
CMolecule::atom_indices(const std::string &label) const -> std::vector<int>
{
    std::vector<int> indices;

    indices.reserve(_identifiers.size());

    const auto identifier = chem_elem::identifier(format::upper_case(label));

    std::ranges::for_each(std::views::iota(0, number_of_atoms()), [&](const int i) {
        if (_identifiers[i] == identifier) indices.push_back(i);
    });

    return indices;
}

auto
CMolecule::nuclear_repulsion_energy() const -> double
{
    double nenergy = 0.0;

    std::ranges::for_each(views::upper_triangular(_identifiers.size()), [&](const auto &index) {
        const auto [i, j] = index;
        nenergy += _identifiers[i] * _identifiers[j] / _coordinates[i].distance(_coordinates[j]);
    });

    return nenergy;
}

auto
CMolecule::check_proximity(const double distance) const -> bool
{
    const auto r2dist = distance * distance;

    auto tru_ij = views::upper_triangular(_identifiers.size());

    auto it = std::ranges::find_if(tru_ij, [&](const auto &index) {
        const auto [i, j] = index;
        return _coordinates[i].distance_square(_coordinates[j]) <= r2dist;
    });

    return it == tru_ij.end();
}

auto
CMolecule::min_distances() const -> std::vector<double>
{
    std::vector<double> mdists(_identifiers.size(), std::numeric_limits<double>::max());

    std::ranges::for_each(views::upper_triangular(_identifiers.size()), [&](const auto &index) {
        const auto [i, j] = index;
        const auto rab    = _coordinates[i].distance(_coordinates[j]);
        if (rab < mdists[i]) mdists[i] = rab;
        if (rab < mdists[j]) mdists[j] = rab;
    });

    return mdists;
}

auto
CMolecule::get_vdw_radii() const -> std::vector<double>
{
    std::vector<double> atomradii;

    auto radii = atomicradii::buildVdwRadii();

    for (int i = 0; i < number_of_atoms(); i++)
    {
        atomradii.push_back(radii[_identifiers[i]]);
    }

    return atomradii;
}

auto
CMolecule::get_mk_radii() const -> std::vector<double>
{
    std::vector<double> atomradii;

    auto mk_radii = atomicradii::buildMkRadii();

    auto vdw_radii = atomicradii::buildVdwRadii();

    for (int i = 0; i < number_of_atoms(); i++)
    {
        if (_identifiers[i] < static_cast<int>(mk_radii.size()))
        {
            atomradii.push_back(mk_radii[_identifiers[i]]);
        }
        else
        {
            atomradii.push_back(vdw_radii[_identifiers[i]]);
        }
    }

    return atomradii;
}

auto
CMolecule::get_chelpg_radii() const -> std::vector<double>
{
    std::vector<double> atomradii;

    auto chelpg_radii = atomicradii::buildChelpgRadii();

    for (int i = 0; i < number_of_atoms(); i++)
    {
        if (_identifiers[i] < static_cast<int>(chelpg_radii.size()))
        {
            atomradii.push_back(chelpg_radii[_identifiers[i]]);
        }
        else
        {
            atomradii.push_back(2.0 / units::bohr_in_angstrom());
        }
    }

    return atomradii;
}

auto
CMolecule::get_covalent_radii() const -> std::vector<double>
{
    std::vector<double> atomradii;

    auto radii = atomicradii::buildCovalentRadii();

    for (int i = 0; i < number_of_atoms(); i++)
    {
        atomradii.push_back(radii[_identifiers[i]]);
    }

    return atomradii;
}

auto
CMolecule::_is_angstrom(const std::string &unit) const -> bool
{
    if (unit.length() >= 3)
    {
        const auto label = std::string("ANGSTROM").substr(0, unit.length());

        return format::upper_case(unit) == label;
    }
    else
    {
        return false;
    }
}
