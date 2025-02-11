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

#ifndef XtbDriver_hpp
#define XtbDriver_hpp

#include <cstdint>
#include <string>
#include <vector>

#ifdef ENABLE_XTB
#include <xtb.h>
#endif

#include "Molecule.hpp"

/**
 Class CXtbDriver enables DFT-B computations using XTB package from Grimme group.
 */
class CXtbDriver
{
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
     Creates a XTB driver object.
     */
    CXtbDriver();

    /**
     Destroys a XTB driver object.
     */
    ~CXtbDriver();

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
    void setMethod(const std::string& method);

    /**
     Sets output filename.

     @param filename the output filename.
     */
    void setOutputFilename(const std::string& filename);

    /**
     Mutes output.
     */
    void mute();

    /**
     Unmutes output.
     */
    void unmute();

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
     Gets XTB method.

     @return the XTB method.
    */
    std::string getMethod() const;

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

    /**
     Gets partial charges as vector (order: natoms).

     @return the partial charges.
    */
    std::vector<double> getPartialCharges() const;

    /**
     Gets bond orders as vector (order: natoms).

     @return the bond orders.
    */
    std::vector<double> getBondOrders() const;

    /**
     Gets number of Atoms.

     @return the number of Atoms.
    */
    int getNumberOfAtoms() const;

    /**
     Gets number of AOs.

     @return the number of AOs.
    */
    int getNumberOfAOs() const;

    /**
     Gets orbital energies as vector (order: naos).

     @return the orbital energies.
    */
    std::vector<double> getOrbitalEnergies() const;

    /**
     Gets orbital occupation numbers as vector (order: naos).

     @return the orbital occupation numbers.
    */
    std::vector<double> getOrbitalOccupations() const;
};

#endif /* XtbDriver_hpp */
