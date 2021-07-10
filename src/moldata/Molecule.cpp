//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include <cmath>
#include <sstream>

#include <mpi.h>

#include "Codata.hpp"
#include "MathFunc.hpp"
#include "StringFormat.hpp"

CMolecule::CMolecule()

    : _charge(0.0)

    , _multiplicity(1)
{
}

CMolecule::CMolecule(const std::vector<double>&      atomCoordinates,
                     const std::vector<double>&      atomCharges,
                     const std::vector<double>&      atomMasses,
                     const std::vector<std::string>& atomLabels,
                     const std::vector<int32_t>&     idsElemental)

    : _charge(0.0)

    , _multiplicity(1)
{
    auto natoms = static_cast<int32_t>(idsElemental.size());

    // set up atom's properties

    _atomCoordinates = CMemBlock2D<double>(atomCoordinates, natoms, 3);

    _atomCharges = CMemBlock<double>(atomCharges);

    _atomMasses = CMemBlock<double>(atomMasses);

    _atomLabels = atomLabels;

    _idsElemental = CMemBlock<int32_t>(idsElemental);

    // set up default indexing of atoms in molecule

    setAtomicIndexes(0);
}

CMolecule::CMolecule(const CMolecule& source)

    : _charge(source._charge)

    , _multiplicity(source._multiplicity)

    , _atomCoordinates(source._atomCoordinates)

    , _atomCharges(source._atomCharges)

    , _atomMasses(source._atomMasses)

    , _atomLabels(source._atomLabels)

    , _idsAtomic(source._idsAtomic)

    , _idsElemental(source._idsElemental)
{
}

CMolecule::CMolecule(CMolecule&& source) noexcept

    : _charge(std::move(source._charge))

    , _multiplicity(std::move(source._multiplicity))

    , _atomCoordinates(std::move(source._atomCoordinates))

    , _atomCharges(std::move(source._atomCharges))

    , _atomMasses(std::move(source._atomMasses))

    , _atomLabels(std::move(source._atomLabels))

    , _idsAtomic(std::move(source._idsAtomic))

    , _idsElemental(std::move(source._idsElemental))
{
}

CMolecule::CMolecule(const CMolecule& mol_1, const CMolecule& mol_2)
{
    std::vector<double>      atomCoordinates;
    std::vector<double>      atomCharges;
    std::vector<double>      atomMasses;
    std::vector<std::string> atomLabels;
    std::vector<int32_t>     idsElemental;

    // x, y, and z coordinates

    for (int32_t i = 0; i < mol_1.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_1._atomCoordinates.data(0)[i]);
    }

    for (int32_t i = 0; i < mol_2.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_2._atomCoordinates.data(0)[i]);
    }

    for (int32_t i = 0; i < mol_1.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_1._atomCoordinates.data(1)[i]);
    }

    for (int32_t i = 0; i < mol_2.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_2._atomCoordinates.data(1)[i]);
    }

    for (int32_t i = 0; i < mol_1.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_1._atomCoordinates.data(2)[i]);
    }

    for (int32_t i = 0; i < mol_2.getNumberOfAtoms(); i++)
    {
        atomCoordinates.push_back(mol_2._atomCoordinates.data(2)[i]);
    }

    // charges, masses, labels, ids

    for (int32_t i = 0; i < mol_1.getNumberOfAtoms(); i++)
    {
        atomCharges.push_back(mol_1._atomCharges.data()[i]);

        atomMasses.push_back(mol_1._atomMasses.data()[i]);

        atomLabels.push_back(mol_1._atomLabels.data()[i]);

        idsElemental.push_back(mol_1._idsElemental.data()[i]);
    }

    for (int32_t i = 0; i < mol_2.getNumberOfAtoms(); i++)
    {
        atomCharges.push_back(mol_2._atomCharges.data()[i]);

        atomMasses.push_back(mol_2._atomMasses.data()[i]);

        atomLabels.push_back(mol_2._atomLabels.data()[i]);

        idsElemental.push_back(mol_2._idsElemental.data()[i]);
    }

    // set up combined molecule

    auto natoms = static_cast<int32_t>(idsElemental.size());

    _atomCoordinates = CMemBlock2D<double>(atomCoordinates, natoms, 3);

    _atomCharges = CMemBlock<double>(atomCharges);

    _atomMasses = CMemBlock<double>(atomMasses);

    _atomLabels = atomLabels;

    _idsElemental = CMemBlock<int32_t>(idsElemental);

    setAtomicIndexes(0);

    setCharge(mol_1._charge + mol_2._charge);

    setMultiplicity(mol_1._multiplicity + mol_2._multiplicity - 1);
}

CMolecule::~CMolecule()
{
}

CMolecule
CMolecule::getSubMolecule(int32_t startIndex, int32_t numAtoms)
{
    std::vector<double>      atomCoordinates;
    std::vector<double>      atomCharges;
    std::vector<double>      atomMasses;
    std::vector<std::string> atomLabels;
    std::vector<int32_t>     idsElemental;

    // boundary check

    auto total_natoms = getNumberOfAtoms();

    if ((startIndex < 0) || (numAtoms <= 0) || (startIndex + numAtoms > total_natoms))
    {
        return CMolecule();
    }

    // x, y, and z coordinates

    for (int32_t i = startIndex; i < startIndex + numAtoms; i++)
    {
        atomCoordinates.push_back(_atomCoordinates.data(0)[i]);
    }

    for (int32_t i = startIndex; i < startIndex + numAtoms; i++)
    {
        atomCoordinates.push_back(_atomCoordinates.data(1)[i]);
    }

    for (int32_t i = startIndex; i < startIndex + numAtoms; i++)
    {
        atomCoordinates.push_back(_atomCoordinates.data(2)[i]);
    }

    // charges, masses, labels, ids

    for (int32_t i = startIndex; i < startIndex + numAtoms; i++)
    {
        atomCharges.push_back(_atomCharges.data()[i]);

        atomMasses.push_back(_atomMasses.data()[i]);

        atomLabels.push_back(_atomLabels.data()[i]);

        idsElemental.push_back(_idsElemental.data()[i]);
    }

    // create sub-molecule

    return CMolecule(atomCoordinates, atomCharges, atomMasses, atomLabels, idsElemental);
}

CMolecule&
CMolecule::operator=(const CMolecule& source)
{
    if (this == &source) return *this;

    _charge = source._charge;

    _multiplicity = source._multiplicity;

    _atomCoordinates = source._atomCoordinates;

    _atomCharges = source._atomCharges;

    _atomMasses = source._atomMasses;

    _atomLabels = source._atomLabels;

    _idsAtomic = source._idsAtomic;

    _idsElemental = source._idsElemental;

    return *this;
}

CMolecule&
CMolecule::operator=(CMolecule&& source) noexcept
{
    if (this == &source) return *this;

    _charge = std::move(source._charge);

    _multiplicity = std::move(source._multiplicity);

    _atomCoordinates = std::move(source._atomCoordinates);

    _atomCharges = std::move(source._atomCharges);

    _atomMasses = std::move(source._atomMasses);

    _atomLabels = std::move(source._atomLabels);

    _idsAtomic = std::move(source._idsAtomic);

    _idsElemental = std::move(source._idsElemental);

    return *this;
}

bool
CMolecule::operator==(const CMolecule& other) const
{
    if (std::fabs(_charge - other._charge) > 1.0e-13) return false;

    if (_multiplicity != other._multiplicity) return false;

    if (_atomCoordinates != other._atomCoordinates) return false;

    if (_atomCharges != other._atomCharges) return false;

    if (_atomMasses != other._atomMasses) return false;

    if (_atomLabels.size() != other._atomLabels.size()) return false;

    for (size_t i = 0; i < _atomLabels.size(); i++)
    {
        if (_atomLabels[i] != other._atomLabels[i]) return false;
    }

    if (_idsAtomic != other._idsAtomic) return false;

    if (_idsElemental != other._idsElemental) return false;

    return true;
}

bool
CMolecule::operator!=(const CMolecule& other) const
{
    return !(*this == other);
}

void
CMolecule::setAtomicIndexes(const int32_t startIndex)
{
    for (int32_t i = 0; i < _idsAtomic.size(); i++)
    {
        _idsAtomic.at(i) = startIndex + i;
    }
}

void
CMolecule::setCharge(const double charge)
{
    _charge = charge;
}

void
CMolecule::setMultiplicity(const int32_t multiplicity)
{
    _multiplicity = multiplicity;
}

double
CMolecule::getCharge() const
{
    return _charge;
}

int32_t
CMolecule::getMultiplicity() const
{
    return _multiplicity;
}

int32_t
CMolecule::getNumberOfAtoms() const
{
    return _idsElemental.size();
}

int32_t
CMolecule::getNumberOfAtoms(const int32_t idElemental) const
{
    return getNumberOfAtoms(0, _idsElemental.size(), idElemental);
}

int32_t
CMolecule::getNumberOfAtoms(const int32_t iAtom, const int32_t nAtoms, const int32_t idElemental) const
{
    int32_t natoms = 0;

    for (int32_t i = iAtom; i < (iAtom + nAtoms); i++)
    {
        if (_idsElemental.at(i) == idElemental) natoms++;
    }

    return natoms;
}

std::set<int32_t>
CMolecule::getElementalComposition() const
{
    std::set<int32_t> elemset;

    for (int32_t i = 0; i < _idsElemental.size(); i++)
    {
        elemset.insert(_idsElemental.at(i));
    }

    return elemset;
}

int32_t
CMolecule::getNumberOfElectrons() const
{
    double nelectrons = -_charge;

    for (int32_t i = 0; i < _atomCharges.size(); i++)
    {
        nelectrons += _atomCharges.at(i);
    }

    return static_cast<int32_t>(nelectrons);
}

const int32_t*
CMolecule::getIdsElemental() const
{
    return _idsElemental.data();
}

const double*
CMolecule::getCoordinatesX() const
{
    return _atomCoordinates.data(0);
}

const double*
CMolecule::getCoordinatesY() const
{
    return _atomCoordinates.data(1);
}

const double*
CMolecule::getCoordinatesZ() const
{
    return _atomCoordinates.data(2);
}

CMemBlock2D<double>
CMolecule::getCoordinates() const
{
    return _atomCoordinates;
}

CMemBlock<double>
CMolecule::getCharges() const
{
    return _atomCharges;
}

CMemBlock<double>
CMolecule::getMasses() const
{
    return _atomMasses;
}

CMemBlock<double>
CMolecule::getMinDistances() const
{
    // allocate and initialize distances

    auto natoms = getNumberOfAtoms();

    CMemBlock<double> mdists(natoms);

    auto rmin = mdists.data();

    mathfunc::set_to(rmin, 1.0e24, natoms);

    // set pointers to coordinates

    auto rx = getCoordinatesX();

    auto ry = getCoordinatesY();

    auto rz = getCoordinatesZ();

    // determine distance to closest neighbouring atom

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto rab = mathfunc::distance(rx[i], ry[i], rz[i], rx[j], ry[j], rz[j]);

            if (rab < rmin[i]) rmin[i] = rab;

            if (rab < rmin[j]) rmin[j] = rab;
        }
    }

    return mdists;
}

double
CMolecule::getNuclearRepulsionEnergy() const
{
    // set up pointers to atomic coordinates anf charges

    auto coordx = getCoordinatesX();

    auto coordy = getCoordinatesY();

    auto coordz = getCoordinatesZ();

    // loop over atoms

    double enuc = 0;

    auto natoms = getNumberOfAtoms();

    for (int32_t i = 0; i < natoms; i++)
    {
        auto rax = coordx[i];

        auto ray = coordy[i];

        auto raz = coordz[i];

        auto zea = _atomCharges.at(i);

        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto rab = mathfunc::distance(rax, ray, raz, coordx[j], coordy[j], coordz[j]);

            enuc += zea * _atomCharges.at(j) / rab;
        }
    }

    return enuc;
}

std::vector<double>
CMolecule::getVdwRadii() const
{
    std::vector<double> atomradii;

    auto radii = vdwradii::buildVdwRadii();

    for (int32_t i = 0; i < getNumberOfAtoms(); i++)
    {
        atomradii.push_back(radii[_idsElemental.data()[i]]);
    }

    return atomradii;
}

std::string
CMolecule::getLabel(const int32_t iAtom) const
{
    if (iAtom < getNumberOfAtoms())
    {
        return _atomLabels[iAtom];
    }

    return std::string();
}

std::string
CMolecule::printGeometry() const
{
    std::stringstream ss;

    ss << "Molecular Geometry (Angstroms)\n";

    ss << std::string(32, '=') << "\n\n";

    ss << "  Atom ";

    ss << fstr::format(std::string("Coordinate X"), 20, fmt::right);

    ss << "  ";

    ss << fstr::format(std::string("Coordinate Y"), 20, fmt::right);

    ss << "  ";

    ss << fstr::format(std::string("Coordinate Z"), 20, fmt::right);

    ss << "  \n\n";

    auto factor = units::getBohrValueInAngstroms();

    auto coordx = _atomCoordinates.data(0);

    auto coordy = _atomCoordinates.data(1);

    auto coordz = _atomCoordinates.data(2);

    for (int32_t i = 0; i < _atomCoordinates.size(0); i++)
    {
        std::string label("  ");

        label.append(_atomLabels.at(i));

        ss << fstr::format(label, 6, fmt::left);

        ss << fstr::to_string(factor * coordx[i], 12, 22, fmt::right);

        ss << fstr::to_string(factor * coordy[i], 12, 22, fmt::right);

        ss << fstr::to_string(factor * coordz[i], 12, 22, fmt::right);

        ss << "\n";
    }

    ss << "\n";

    return ss.str();
}

bool
CMolecule::checkProximity(const double minDistance) const
{
    auto natoms = _atomCoordinates.size(0);

    auto coordx = _atomCoordinates.data(0);

    auto coordy = _atomCoordinates.data(1);

    auto coordz = _atomCoordinates.data(2);

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto rab = mathfunc::distance(coordx[i], coordy[i], coordz[i], coordx[j], coordy[j], coordz[j]);

            if (rab < minDistance)
            {
                return false;
            }
        }
    }

    return true;
}

void
CMolecule::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_charge, comm);

        mpi::bcast(_multiplicity, comm);

        _atomCoordinates.broadcast(rank, comm);

        _atomCharges.broadcast(rank, comm);

        _atomMasses.broadcast(rank, comm);

        mpi::bcast(_atomLabels, rank, comm);

        _idsAtomic.broadcast(rank, comm);

        _idsElemental.broadcast(rank, comm);
    }
}

std::ostream&
operator<<(std::ostream& output, const CMolecule& source)
{
    output << std::endl;

    output << "[CMolecule (Object):" << &source << "]" << std::endl;

    output << "_charge: " << source._charge << std::endl;

    output << "_multiplicity: " << source._multiplicity << std::endl;

    output << "_atomCoordinates: " << source._atomCoordinates << std::endl;

    output << "_atomCharges: " << source._atomCharges << std::endl;

    output << "_atomMasses: " << source._atomMasses << std::endl;

    output << "_atomLabels: " << std::endl;

    for (size_t i = 0; i < source._atomLabels.size(); i++)
    {
        output << "atomsLabels_[" << i << "]" << std::endl;

        output << source._atomLabels[i] << std::endl;
    }

    output << "_idsAtomic: " << source._idsAtomic << std::endl;

    output << "_idsElemental: " << source._idsElemental << std::endl;

    return output;
}
