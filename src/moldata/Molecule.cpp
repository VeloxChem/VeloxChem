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
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "MathFunc.hpp"
#include "StringFormat.hpp"

CMolecule::CMolecule(const std::vector<int64_t>& identifiers, const std::vector<TPoint3D>& coordinates, const std::string& units)
{
    if (const auto natoms = identifiers.size(); coordinates.size() == natoms)
    {
        for (size_t i = 0; i < natoms; i++)
        {
            addAtom(identifiers[i], coordinates[i], units);
        }
    }
}

CMolecule::CMolecule(const std::vector<std::string>& labels, const std::vector<TPoint3D>& coordinates, const std::string& units)
{
    if (const auto natoms = labels.size(); coordinates.size() == natoms)
    {
        for (size_t i = 0; i < natoms; i++)
        {
            addAtom(labels[i], coordinates[i], units);
        }
    }
}

CMolecule::CMolecule(const CMolecule& molfrag_one, const CMolecule& molfrag_two)
{
    // copy geometrical data from molecule A

    _identifiers = molfrag_one._identifiers;

    _coordinates = molfrag_one._coordinates;

    // append geometrical data from molecule B

    if (const auto natoms = molfrag_two.getNumberOfAtoms(); natoms > 0)
    {
        for (int64_t i = 0; i < natoms; i++)
        {
            _identifiers.push_back(molfrag_two._identifiers[i]);

            _coordinates.push_back(molfrag_two._coordinates[i]);
        }
    }

    // set charge and spin state

    _charge = molfrag_one._charge + molfrag_two._charge;

    const auto spin_one = (molfrag_one._multiplicity - 1) / 2;

    const auto spin_two = (molfrag_two._multiplicity - 1) / 2;

    _multiplicity = 2 * (spin_one + spin_two) + 1;
}

auto
CMolecule::addAtom(const std::string& label, const TPoint3D& coordinates, const std::string& units) -> void
{
    if (CChemicalElement elem; elem.setAtomType(fstr::upcase(label)))
    {
        addAtom(elem.getIdentifier(), coordinates, units);
    }
}

auto
CMolecule::addAtom(const int64_t identifier, const TPoint3D& coordinates, const std::string& units) -> void
{
    if (CChemicalElement elem; (identifier >= 0) && (identifier <= elem.getMaxIdentifier()))
    {
        _identifiers.push_back(identifier);

        if (_isAngstroms(units))
        {
            const auto fact = 1.0 / units::getBohrValueInAngstroms();

            _coordinates.push_back({coordinates[0] * fact, coordinates[1] * fact, coordinates[2] * fact});
        }
        else
        {
            _coordinates.push_back(coordinates);
        }
    }
    else
    {
        std::cerr << "*** Unsupported chemical element with elemental number ";

        std::cerr << identifier << " is encountered!!!" << std::endl;

        std::exit(EXIT_FAILURE);
    }
}

auto
CMolecule::setCharge(const double charge) -> void
{
    _charge = charge;
}

auto
CMolecule::setMultiplicity(const int64_t multiplicity) -> void
{
    _multiplicity = multiplicity;
}

auto
CMolecule::getCharge() const -> double
{
    return _charge;
}

auto
CMolecule::getMultiplicity() const -> int64_t
{
    return _multiplicity;
}

auto
CMolecule::getNumberOfAtoms() const -> int64_t
{
    return static_cast<int64_t>(_identifiers.size());
}

auto
CMolecule::getNumberOfAtoms(const int64_t identifier) const -> int64_t
{
    int64_t count = 0;

    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        for (int64_t i = 0; i < natoms; i++)
        {
            if (_identifiers[i] == identifier) count++;
        }
    }

    return count;
}

auto
CMolecule::getNumberOfAtoms(const int64_t iatom, const int64_t natoms, const int64_t identifier) const -> int64_t
{
    int64_t count = 0;

    if (natoms > 0)
    {
        for (int64_t i = iatom; i < (iatom + natoms); i++)
        {
            if (_identifiers[i] == identifier) count++;
        }
    }

    return count;
}

auto
CMolecule::getElementalComposition() const -> std::set<int64_t>
{
    std::set<int64_t> elem_comp;

    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        for (int64_t i = 0; i < natoms; i++)
        {
            elem_comp.insert(_identifiers[i]);
        }
    }

    return elem_comp;
}

auto
CMolecule::getNumberOfElectrons() const -> int64_t
{
    auto qsum = -_charge;

    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        auto charges = getCharges();

        for (int64_t i = 0; i < natoms; i++)
        {
            qsum += charges[i];
        }
    }

    return qsum;
}

auto
CMolecule::getIdsElemental() const -> std::vector<int64_t>
{
    return _identifiers;
}

auto
CMolecule::getCoordinates(const std::string& units) const -> std::vector<TPoint3D>
{
    if (_isAngstroms(units))
    {
        if (const auto natoms = getNumberOfAtoms(); natoms > 0)
        {
            std::vector<TPoint3D> coords;

            const auto fact = units::getBohrValueInAngstroms();

            for (int64_t i = 0; i < natoms; i++)
            {
                const auto rxyz = _coordinates[i];

                coords.push_back({rxyz[0] * fact, rxyz[1] * fact, rxyz[2] * fact});
            }

            return coords;
        }
        else
        {
            return _coordinates;
        }
    }
    else
    {
        return _coordinates;
    }
}

auto
CMolecule::getCharges() const -> std::vector<double>
{
    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        std::vector<double> charges;

        for (int64_t i = 0; i < natoms; i++)
        {
            if (CChemicalElement elem; elem.setAtomType(_identifiers[i]))
            {
                charges.push_back(elem.getAtomicCharge());
            }
            else
            {
                charges.push_back(0.0);
            }
        }

        return charges;
    }
    else
    {
        return std::vector<double>();
    }
}

auto
CMolecule::getMasses() const -> std::vector<double>
{
    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        std::vector<double> masses;

        for (int64_t i = 0; i < natoms; i++)
        {
            if (CChemicalElement elem; elem.setAtomType(_identifiers[i]))
            {
                masses.push_back(elem.getAtomicMass());
            }
            else
            {
                masses.push_back(0.0);
            }
        }

        return masses;
    }
    else
    {
        return std::vector<double>();
    }
}

auto
CMolecule::getLabels() const -> std::vector<std::string>
{
    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        std::vector<std::string> labels;

        for (int64_t i = 0; i < natoms; i++)
        {
            labels.push_back(getLabel(i));
        }

        return labels;
    }
    else
    {
        return std::vector<std::string>();
    }
}

auto
CMolecule::getLabel(const int64_t iatom) const -> std::string
{
    if (CChemicalElement elem; elem.setAtomType(_identifiers[iatom]))
    {
        return elem.getName();
    }
    else
    {
        return std::string();
    }
}

auto
CMolecule::getAtomCoordinates(const int64_t iatom, const std::string& units) const -> TPoint3D
{
    if (_isAngstroms(units))
    {
        const auto fact = units::getBohrValueInAngstroms();

        const auto rxyz = _coordinates[iatom];

        return {rxyz[0] * fact, rxyz[1] * fact, rxyz[2] * fact};
    }
    else
    {
        return _coordinates[iatom];
    }
}

auto
CMolecule::getAtomIndexes(const std::string& label) const -> std::vector<int64_t>
{
    if (CChemicalElement elem; elem.setAtomType(label))
    {
        if (const auto natoms = getNumberOfAtoms(); natoms > 0)
        {
            std::vector<int64_t> indexes;

            const auto idatm = elem.getIdentifier();

            for (int64_t i = 0; i < natoms; i++)
            {
                if (idatm == _identifiers[i]) indexes.push_back(i);
            }

            return indexes;
        }
        else
        {
            return std::vector<int64_t>();
        }
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecule::getNuclearRepulsionEnergy() const -> double
{
    double nenergy = 0.0;

    if (const auto natoms = getNumberOfAtoms(); natoms > 1)
    {
        const auto charges = getCharges();

        for (int64_t i = 0; i < natoms; i++)
        {
            if (_identifiers[i] != 0)
            {
                const auto [ax, ay, az] = _coordinates[i];

                const auto zea = charges[i];

                for (int64_t j = i + 1; j < natoms; j++)
                {
                    if (_identifiers[j] != 0)
                    {
                        const auto [bx, by, bz] = _coordinates[j];

                        const auto rabx = ax - bx;

                        const auto raby = ay - by;

                        const auto rabz = az - bz;

                        nenergy += zea * charges[j] / std::sqrt(rabx * rabx + raby * raby + rabz * rabz);
                    }
                }
            }
        }
    }

    return nenergy;
}

auto
CMolecule::checkProximity(const double minDistance) const -> bool
{
    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        const auto r2min = minDistance * minDistance;

        for (int64_t i = 0; i < natoms; i++)
        {
            if (_identifiers[i] != 0)
            {
                const auto [ax, ay, az] = _coordinates[i];

                for (int64_t j = i + 1; j < natoms; j++)
                {
                    if (_identifiers[j] != 0)
                    {
                        const auto [bx, by, bz] = _coordinates[j];

                        const auto rabx = ax - bx;

                        const auto raby = ay - by;

                        const auto rabz = az - bz;

                        const auto r2ab = rabx * rabx + raby * raby + rabz * rabz;

                        if (r2ab < r2min) return false;
                    }
                }
            }
        }
    }

    return true;
}

auto
CMolecule::printGeometry() const -> std::string
{
    std::stringstream ss;

    ss << "Molecular Geometry (Angstroms)\n";

    ss << std::string(32, '=') << "\n\n";

    ss << "  Atom ";

    ss << fstr::format(std::string("Coordinate X"), 20, fmt_t::right);

    ss << "  ";

    ss << fstr::format(std::string("Coordinate Y"), 20, fmt_t::right);

    ss << "  ";

    ss << fstr::format(std::string("Coordinate Z"), 20, fmt_t::right);

    ss << "  \n\n";

    if (const auto natoms = getNumberOfAtoms(); natoms > 0)
    {
        const auto labels = getLabels();

        const auto coords = getCoordinates("angstrom");

        for (int64_t i = 0; i < natoms; i++)
        {
            if (_identifiers[i] != 0)
            {
                std::string label("  ");

                label.append(labels[i]);

                ss << fstr::format(label, 6, fmt_t::left);

                const auto [rx, ry, rz] = coords[i];

                ss << fstr::to_string(rx, 12, 22, fmt_t::right);

                ss << fstr::to_string(ry, 12, 22, fmt_t::right);

                ss << fstr::to_string(rz, 12, 22, fmt_t::right);

                ss << "\n";
            }
        }

        ss << "\n";
    }

    return ss.str();
}

auto
CMolecule::_isAngstroms(const std::string& units) const -> bool
{
    return (units.length() >= 3) && (fstr::upcase(units) == std::string("ANGSTROM").substr(0, units.length()));
}

auto
CMolecule::getMinDistances() const -> std::vector<double>
{
    // allocate and initialize distances

    auto natoms = getNumberOfAtoms();

    std::vector<double> mdists(natoms);

    auto rmin = mdists.data();

    std::fill(rmin, rmin + natoms, 1.0e24);

    // set pointers to coordinates

    auto coords = getCoordinates(std::string("BOHR"));

    // determine distance to closest neighbouring atom

    for (int64_t i = 0; i < natoms; i++)
    {
        const auto& ri = coords[i];

        for (int64_t j = i + 1; j < natoms; j++)
        {
            const auto& rj = coords[j];

            auto rab = mathfunc::distance(ri[0], ri[1], ri[2], rj[0], rj[1], rj[2]);

            if (rab < rmin[i]) rmin[i] = rab;

            if (rab < rmin[j]) rmin[j] = rab;
        }
    }

    return mdists;
}
