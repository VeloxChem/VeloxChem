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

#ifndef XTBDriver_hpp
#define XTBDriver_hpp

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

#ifdef ENABLE_XTB
#include <xtb.h>
#endif

class CMolecule;

/**
 Class CXTBDriver enables DFT-B computations using XTB package from Grimme group.

 @author Z. Rinkevicius
 */
class CXTBDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     The name of the XTB output file.
     */
    std::string _outputFilename;

    /**
     The XTB method.
     */
    std::string _xtbMethod;

    /**
    The electronic temperature.
    */
    double _electronicTemp;

    /**
    The maximum number of SCF iterations.
    */
    int _maxIterations;

    /**
     The number of atoms in XTB molecular system.
    */
    int _natoms;

#ifdef ENABLE_XTB
    /**
     The XTB package environment.
    */
    xtb_TEnvironment _environment;

    /**
     The DFT-B calculator.
    */
    xtb_TCalculator _calculator;

    /**
     The DFT-B results data.
    */
    xtb_TResults _results;

    /**
     Converts VeloxChem molecule to XTB molecular system data structure.

     @param molecule the molecule.
     @returm the XTB  molecular system data structure.
    */
    xtb_TMolecule _set_molecule(const CMolecule& molecule);
#endif

   public:
    /**
     Creates a XTB driver object using MPI info.

     @param comm the MPI communicator.
     */
    CXTBDriver(MPI_Comm comm);

    /**
     Destroys a XTB driver object.
     */
    ~CXTBDriver();

    /**
     Sets maximum number of SCF iterations.

     @param maxIterations the maximum number of iterations.
     */
    void setMaxIterations(const int maxIterations);

    /**
     Sets electronic temperature for electron smearing.

     @param electronicTemp the electronic temperature in Kelvins.
     */
    void setElectronicTemp(const double electronicTemp);

    /**
     Sets XTB method.

     @param method the XTB method.
     */
    void setMethod(const std::string method);

    /**
     Sets output filename.

     @param filename the output filename.
     */
    void setOutputFilename(const std::string filename);

    /**
     Computes DTB-B single point energy using XTB package.

     @param molecule the molecule.
     */
    void compute(const CMolecule& molecule);

    /**
     Checks if XTB package is available.

     @return true if XTB package available, false otherwise.
     */
    static constexpr auto
    isAvailable() -> bool
    {
#ifdef ENABLE_XTB
        return true;
#endif
        return false;
    }

    /**
     Checks if XTB driver is running on master node.

     @return true if XTB driver is running on master node, false otherwise.
     */
    bool isMasterNode() const;

    /**
     Gets state of XTB environment.

     @return the state of XTB environment.
    */
    bool getState() const;

    /**
     Gets XTB output as a vector of strings.

     @return a vector of strings.
    */
    std::vector<std::string> getOutput() const;

    /**
     Gets XTB output filename.

     @return the output filename.
    */
    std::string getOutputFilename() const;

    /**
     Gets energy of molecular system.

     @return the energy of molecular system.
    */
    double getEnergy() const;

    /**
     Gets molecular gradient as vector (order: natoms x 3).

     @return the molecular gradient.
    */
    std::vector<double> getGradient() const;

    /**
     Gets molecular dipole moment as vector (order: 3).

     @return the molecular dipole moment.
    */
    std::vector<double> getDipole() const;
};

#endif /* XTBDriver_hpp */
