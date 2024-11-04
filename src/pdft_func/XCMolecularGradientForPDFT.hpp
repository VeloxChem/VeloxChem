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

#ifndef XCMolecularGradientForPDFT_hpp
#define XCMolecularGradientForPDFT_hpp

#include "XCPairDensityFunctional.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"

namespace xcgradpdft {  // xcgradpdft namespace
    /**
     Integrates first-order LDA exchange-correlation functional contribution to
     PDFT molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param twoBodyDensityMatrix the MO two-body active density matrix.
     @param activeMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc exchange-correlation functional.
     @param rs_omega range-separation parameter.
     @return the molecular gradient.
     */
     CDenseMatrix integrateVxcPDFTGradientForLDA(const CMolecule&                molecule,
                                                 const CMolecularBasis&          basis,
                                                 const double*                   densityMatrixPointer,
                                                 const CDenseMatrix&             twoBodyDensityMatrix,
                                                 const CDenseMatrix&             activeMOs,
                                                 const CMolecularGrid&           molecularGrid,
                                                 const double                    screeningThresholdForGTOValues,
                                                 const CXCPairDensityFunctional& xcFunctional,
                                                 const double                    rs_omega);
    /**
     Integrates first-order GGA exchange-correlation functional contribution to
     PDFT molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param twoBodyDensityMatrix the MO two-body active density matrix.
     @param activeMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc exchange-correlation functional.
     @param rs_omega range-separation parameter.
     @return the molecular gradient.
     */
     CDenseMatrix integrateVxcPDFTGradientForGGA(const CMolecule&                molecule,
                                                 const CMolecularBasis&          basis,
                                                 const double*                   densityMatrixPointer,
                                                 const CDenseMatrix&             twoBodyDensityMatrix,
                                                 const CDenseMatrix&             activeMOs,
                                                 const CMolecularGrid&           molecularGrid,
                                                 const double                    screeningThresholdForGTOValues,
                                                 const CXCPairDensityFunctional& xcFunctional,
                                                 const double                    rs_omega);

    // Duplicate for now
    void _computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis);
}   // namespace xcgradpdft

#endif /* XCMolecularGradientForPDFT_hpp */
