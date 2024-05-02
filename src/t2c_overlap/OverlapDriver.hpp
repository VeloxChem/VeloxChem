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

#ifndef OverlapDriver_hpp
#define OverlapDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class COverlapDriver provides methods for computing two-center
 overlap integrals.

 @author Z. Rinkevicius
 */
class COverlapDriver
{
   public:
    /**
     Creates an overlap integrals driver.
     */
    COverlapDriver() = default;

    /**
     Computes overlap matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix;
    
    /**
     Computes overlap matrix for given molecule and pair of molecular bases.

     @param bra_basis the molecular basis on bra side.
     @param ket_basis the molecular basis on ket side.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto compute(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis, const CMolecule& molecule) const -> CMatrix;
};

#endif /* OverlapDriver_hpp */
