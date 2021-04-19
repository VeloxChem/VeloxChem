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

#ifndef GtoTransform_hpp
#define GtoTransform_hpp

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtotra {  // gtotra namespace
    
    /**
     Transforms AO matrix from Dalton to VeloxChem format.

     @param matrix the AO matrix in Dalton format.
     @param basis the molecular basis set.
     @param molecule the molecule.
     @return the AO matrix in VeloxChem format.
     */
    CDenseMatrix to_veloxchem(const CDenseMatrix&    matrix,
                              const CMolecularBasis& basis,
                              const CMolecule&       molecule);
    
    /**
     Transforms AO matrix from VeloxChem to Dalton format.
     
     @param matrix the AO matrix in VeloxChem format.
     @param basis the molecular basis set.
     @param molecule the molecule.
     @return the AO matrix in Dalton format.
     */
    CDenseMatrix to_dalton(const CDenseMatrix&    matrix,
                           const CMolecularBasis& basis,
                           const CMolecule&       molecule);
}  // namespace gtotra

#endif /* GtoTransform_hpp */
