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

#ifndef NewDispersionModel_hpp
#define NewDispersionModel_hpp

#include <string>

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

/**
 Class CNewDispersionModel computes dispersion energy and gradient using dftd4.
 */
class CNewDispersionModel
{
    /**
     The dispersion energy.
     */
    double _energy;

    /**
     The dispersion gradient (Nx3).
     */
    CDenseMatrix _gradient;

   public:
    /**
     Creates a dispersion model object.
     */
    CNewDispersionModel();

    /**
     Destroys a dispersion model object.
     */
    ~CNewDispersionModel();

    /**
     Computes dispersion energy and gradient for a given molecule and a given
     density functional.

     @param molecule the molecule.
     @param xcLabel the label of the density functional.
     */
    void compute(const CMolecule& molecule, const std::string& xcLabel);

    void check_error_code(const int error_code, const std::string& msg);

    /**
     Gets dispersion energy.

     @return the dispersion energy.
     */
    double getEnergy() const;

    /**
     Gets dispersion gradient.

     @return the dispersion gradient.
     */
    CDenseMatrix getGradient() const;
};

#endif /* NewDispersionModel_hpp */
