//
//                           VELOXCHEM 1.0-RC
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

#include "XTBDriver.hpp"

#include <iostream>
#include <fstream>

#include <mpi.h>

#include "Molecule.hpp"
#include "MpiFunc.hpp"

CXTBDriver::CXTBDriver(MPI_Comm comm)
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
#endif
}

CXTBDriver::~CXTBDriver()
{
#ifdef ENABLE_XTB
    xtb_delResults(&_results);
    
    xtb_delCalculator(&_calculator);
    
    xtb_delEnvironment(&_environment);
#endif
}

void 
CXTBDriver::setMaxIterations(const int maxIterations)
{
    _maxIterations = maxIterations; 
}

void 
CXTBDriver::setElectronicTemp(const double electronicTemp)
{
    _electronicTemp = electronicTemp; 
}

void 
CXTBDriver::setMethod(const std::string method)
{
    _xtbMethod = method; 
}

void 
CXTBDriver::setOutputFilename(const std::string filename)
{
    _outputFilename = filename; 
}

void 
CXTBDriver::compute(const CMolecule& molecule)
{
#ifdef ENABLE_XTB
    if (isMasterNode())
    {
        // set up output stream 
    
        xtb_setOutput(_environment, _outputFilename.c_str());

        xtb_setVerbosity(_environment, XTB_VERBOSITY_FULL); 

        // xtb_setAccuracy(_environment, _calculator, 0.0001); 

        // set up molecular data structure
    
        auto tmol = _set_molecule(molecule);

        // load DFT-B method parameters
    
        if (_xtbMethod == "gfn2")
        {
            xtb_loadGFN2xTB(_environment, tmol, _calculator, NULL);
        }
        else if (_xtbMethod == "gfn1")
        {
            xtb_loadGFN1xTB(_environment, tmol, _calculator, NULL);
        }
        else if (_xtbMethod == "gfn0")
        {
            xtb_loadGFN0xTB(_environment, tmol, _calculator, NULL);
        }
        else
        {
            return;
        }
  
        // perform single point calculation

        xtb_setMaxIter(_environment, _calculator, _maxIterations);
        
        xtb_setElectronicTemp(_environment, _calculator, _electronicTemp);
        
        xtb_singlepoint(_environment, tmol, _calculator, _results);

        // delete molecular data structure

        xtb_delMolecule(&tmol);

        // release output file
    
        xtb_releaseOutput(_environment);    
    }
#endif
}

bool 
CXTBDriver::isMasterNode() const 
{
    return _locRank == mpi::master(); 
}

bool 
CXTBDriver::isAvailable() const 
{
#ifdef ENABLE_XTB
    return true;
#endif
    return false; 
}

bool
CXTBDriver::getState()
{
#ifdef ENABLE_XTB
    return xtb_checkEnvironment(_environment) > 0;
#endif
    return false;
}

std::vector<std::string>
CXTBDriver::getOutput()
{
    std::vector<std::string> output_strings;

    std::ifstream fst(_outputFilename.c_str());

    if (fst.is_open())
    {
        while(true)
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
CXTBDriver::getOutputFilename()
{
    return _outputFilename;
}

double 
CXTBDriver::getEnergy() 
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
CXTBDriver::getGradient()
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

#ifdef ENABLE_XTB
xtb_TMolecule
CXTBDriver::_set_molecule(const CMolecule& molecule)
{
    _natoms = static_cast<int>(molecule.getNumberOfAtoms());
    
    double charge = molecule.getCharge();
   
    int mult = static_cast<int>(molecule.getMultiplicity());  

    int uhf = (mult > 1) ? mult - 1 : 0;
    
    std::vector<int> atoms(_natoms, 0);
    
    std::vector<double> coords(3 * _natoms, 0.0);
    
    // reformat molecular data
    
    auto eleids = molecule.getIdsElemental();
    
    auto rx =  molecule.getCoordinatesX();
    
    auto ry =  molecule.getCoordinatesY();
    
    auto rz =  molecule.getCoordinatesZ();
    
    for (int i = 0; i < _natoms; i++)
    {
        atoms.at(i) = static_cast<int>(eleids[i]);
        
        coords.at(3 * i) = rx[i];
        
        coords.at(3 * i + 1) = ry[i];
        
        coords.at(3 * i + 2) = rz[i];
    }
    
    return xtb_newMolecule(_environment, &_natoms, atoms.data(),
                           coords.data(),
                           &charge, &uhf,
                           NULL, NULL);
}
#endif
