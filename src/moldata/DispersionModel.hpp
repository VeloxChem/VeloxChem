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

#ifndef DispersionModel_hpp
#define DispersionModel_hpp

#include <string>

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

/**
 Class CDispersionModel computes dispersion energy and gradient using dftd4.
 */
class CDispersionModel
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
    CDispersionModel();

    /**
     Destroys a dispersion model object.
     */
    ~CDispersionModel();

    /**
     Checks if dftd4 is available.

     @return true if dftd4 is available, false otherwise.
     */
    static constexpr bool is_available()
    {
#ifdef ENABLE_DFTD4
        return true;
#else
        return false;
#endif
    }

    /**
     Computes dispersion energy and gradient for a given molecule and a given
     density functional.

     @param molecule the molecule.
     @param xcLabel the label of the density functional.
     */
    void compute(const CMolecule& molecule, const std::string& xcLabel);

    /**
     Checks error code.

     @param error_code the error code.
     @param msg the error message.
     */
    void check_error_code(const int error_code, const std::string& msg) const;

    /**
     Gets dispersion energy.

     @return the dispersion energy.
     */
    double get_energy() const;

    /**
     Gets dispersion gradient.

     @return the dispersion gradient.
     */
    CDenseMatrix get_gradient() const;
};

#endif /* DispersionModel_hpp */
