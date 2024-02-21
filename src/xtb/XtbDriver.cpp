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

#include "XtbDriver.hpp"

#include <mpi.h>

#include <fstream>
#include <iostream>

#include "Molecule.hpp"
#include "MpiFunc.hpp"

CXtbDriver::CXtbDriver(MPI_Comm comm)
    : _outputFilename(std::string("STDOUT"))

    , _xtbMethod(std::string("gfn2"))
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;

    _electronicTemp = 300.0;

    _maxIterations = 280;

    _natoms = 0;

#ifdef ENABLE_XTB
    _environment = xtb_newEnvironment();

    _calculator = xtb_newCalculator();

    _results = xtb_newResults();

    xtb_setVerbosity(_environment, XTB_VERBOSITY_FULL);
#endif
}

CXtbDriver::~CXtbDriver()
{
#ifdef ENABLE_XTB
    xtb_delResults(&_results);

    xtb_delCalculator(&_calculator);

    xtb_delEnvironment(&_environment);
#endif
}

void
CXtbDriver::setMaxIterations(const int maxIterations)
{
    _maxIterations = maxIterations;
}

void
CXtbDriver::setElectronicTemp(const double electronicTemp)
{
    _electronicTemp = electronicTemp;
}

void
CXtbDriver::setMethod(const std::string& method)
{
    _xtbMethod = method;
}

void
CXtbDriver::setOutputFilename(const std::string& filename)
{
    _outputFilename = filename;
}

void
CXtbDriver::mute()
{
#ifdef ENABLE_XTB
    if (isMasterNode()) xtb_setVerbosity(_environment, XTB_VERBOSITY_MUTED);
#endif
}

void
CXtbDriver::unmute()
{
#ifdef ENABLE_XTB
    if (isMasterNode()) xtb_setVerbosity(_environment, XTB_VERBOSITY_FULL);
#endif
}

void
CXtbDriver::compute(const CMolecule& molecule)
{
#ifdef ENABLE_XTB
    if (isMasterNode())
    {
        // set up output stream

        xtb_setOutput(_environment, _outputFilename.c_str());

        // xtb_setAccuracy(_environment, _calculator, 0.0001);

        // set up molecular data structure

        auto tmol = _set_molecule(molecule);

        // load DFT-B method parameters

        if (_xtbMethod == "gfn2")
        {
            xtb_loadGFN2xTB(_environment, tmol, _calculator, nullptr);
        }
        else if (_xtbMethod == "gfn1")
        {
            xtb_loadGFN1xTB(_environment, tmol, _calculator, nullptr);
        }
        else if (_xtbMethod == "gfn0")
        {
            xtb_loadGFN0xTB(_environment, tmol, _calculator, nullptr);
        }
        else
        {
            return;
        }

        // perform single point calculation

        xtb_setMaxIter(_environment, _calculator, _maxIterations);

        xtb_setElectronicTemp(_environment, _calculator, _electronicTemp);

        xtb_singlepoint(_environment, tmol, _calculator, _results);

        errors::assertMsgCritical(xtb_checkEnvironment(_environment) == 0, std::string("XtbDriver: Error in XTB calculation"));

        // delete molecular data structure

        xtb_delMolecule(&tmol);

        // release output file

        xtb_releaseOutput(_environment);
    }
#endif
}

bool
CXtbDriver::isMasterNode() const
{
    return _locRank == mpi::master();
}

bool
CXtbDriver::getState() const
{
#ifdef ENABLE_XTB
    return xtb_checkEnvironment(_environment) == 0;
#endif
    return false;
}

std::vector<std::string>
CXtbDriver::getOutput() const
{
    std::vector<std::string> output_strings;

    std::ifstream fst(_outputFilename.c_str());

    if (fst.is_open())
    {
        while (true)
        {
            std::string str;

            std::getline(fst, str);

            if (fst.eof()) break;

            output_strings.push_back(str);
        }

        fst.close();
    }

    return output_strings;
}

std::string
CXtbDriver::getMethod() const
{
    return _xtbMethod;
}

std::string
CXtbDriver::getOutputFilename() const
{
    return _outputFilename;
}

double
CXtbDriver::getEnergy() const
{
    double energy = 0.0;

    if (isMasterNode())
    {
#ifdef ENABLE_XTB
        xtb_getEnergy(_environment, _results, &energy);
#endif
    }

    return energy;
}

std::vector<double>
CXtbDriver::getGradient() const
{
    std::vector<double> grad;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        grad = std::vector<double>(_natoms * 3, 0.0);

        xtb_getGradient(_environment, _results, grad.data());
#endif
    }

    return grad;
}

std::vector<double>
CXtbDriver::getDipole() const
{
    std::vector<double> dipole;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        dipole = std::vector<double>(3, 0.0);

        xtb_getDipole(_environment, _results, dipole.data());
#endif
    }

    return dipole;
}

std::vector<double>
CXtbDriver::getPartialCharges() const
{
    std::vector<double> partial_charges;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        partial_charges = std::vector<double>(_natoms, 0.0);

        xtb_getCharges(_environment, _results, partial_charges.data());
#endif
    }

    return partial_charges;
}

std::vector<double>
CXtbDriver::getBondOrders() const
{
    std::vector<double> bond_orders;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        bond_orders = std::vector<double>(_natoms * _natoms, 0.0);

        xtb_getBondOrders(_environment, _results, bond_orders.data());
#endif
    }

    return bond_orders;
}

int32_t
CXtbDriver::getNumberOfAtoms() const
{
    return _natoms;
}

int32_t
CXtbDriver::getNumberOfAOs() const
{
    int nao = 0;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        xtb_getNao(_environment, _results, &nao);
#endif
    }

    return static_cast<int32_t>(nao);
}

std::vector<double>
CXtbDriver::getOrbitalEnergies() const
{
    std::vector<double> orbital_energies;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        auto nao = getNumberOfAOs();

        orbital_energies = std::vector<double>(nao, 0.0);

        xtb_getOrbitalEigenvalues(_environment, _results, orbital_energies.data());
#endif
    }

    return orbital_energies;
}

std::vector<double>
CXtbDriver::getOrbitalOccupations() const
{
    std::vector<double> orbital_occupations;

    if ((_natoms > 0) && isMasterNode())
    {
#ifdef ENABLE_XTB
        auto nao = getNumberOfAOs();

        orbital_occupations = std::vector<double>(nao, 0.0);

        xtb_getOrbitalOccupations(_environment, _results, orbital_occupations.data());
#endif
    }

    return orbital_occupations;
}

#ifdef ENABLE_XTB
xtb_TMolecule
CXtbDriver::_set_molecule(const CMolecule& molecule)
{
    _natoms = static_cast<int>(molecule.getNumberOfAtoms());

    double charge = molecule.getCharge();

    int mult = static_cast<int>(molecule.getMultiplicity());

    int uhf = (mult > 1) ? mult - 1 : 0;

    std::vector<int> atoms(_natoms, 0);

    std::vector<double> coords(3 * _natoms, 0.0);

    // reformat molecular data

    auto eleids = molecule.getIdsElemental();

    auto rx = molecule.getCoordinatesX();

    auto ry = molecule.getCoordinatesY();

    auto rz = molecule.getCoordinatesZ();

    for (int i = 0; i < _natoms; i++)
    {
        atoms.at(i) = static_cast<int>(eleids[i]);

        coords.at(3 * i) = rx[i];

        coords.at(3 * i + 1) = ry[i];

        coords.at(3 * i + 2) = rz[i];
    }

    return xtb_newMolecule(_environment, &_natoms, atoms.data(), coords.data(), &charge, &uhf, nullptr, nullptr);
}
#endif
