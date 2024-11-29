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

#ifndef OldOneElecIntsDrivers_hpp
#define OldOneElecIntsDrivers_hpp

#include "AngularMomentumIntegrals.hpp"
#include "LinearMomentumIntegrals.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

namespace oldonee {  // oldonee namespace

class COldOneElecIntsMatrix
{
    std::vector<CDenseMatrix> _matrices;

   public:
    COldOneElecIntsMatrix(const std::vector<CDenseMatrix> matrices);

    CDenseMatrix get_matrix(const int idx);

    ~COldOneElecIntsMatrix();
};

class COldElectricDipoleIntegralsDriver
{
   public:
    COldElectricDipoleIntegralsDriver();

    ~COldElectricDipoleIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

class COldLinearMomentumIntegralsDriver
{
   public:
    COldLinearMomentumIntegralsDriver();

    ~COldLinearMomentumIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

class COldAngularMomentumIntegralsDriver
{
   public:
    COldAngularMomentumIntegralsDriver();

    ~COldAngularMomentumIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

}  // namespace oldonee

#endif /* OldOneElecIntsDrivers_hpp */
