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

#ifndef XCMolecularGradient_hpp
#define XCMolecularGradient_hpp

#include <array>
#include <list>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

/**
 Class CXCMolecularGradient implements XC molecular gradient.
 */
class CXCMolecularGradient
{
   private:
    /**
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

    /**
     Integrates first-order LDA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed AO density matrix (to be contracted
            with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    auto _integrateVxcGradientForLDA(const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const std::vector<const double*>& rwDensityPointers,
                                     const std::vector<const double*>& gsDensityPointers,
                                     const CMolecularGrid&   molecularGrid,
                                     const CXCFunctional&    xcFunctional) const -> CDenseMatrix;

    auto _integrateVxcGradientForLDAOpenShell(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const std::vector<const double*>& rwDensityPointers,
                                              const std::vector<const double*>& gsDensityPointers,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCFunctional&    xcFunctional) const -> CDenseMatrix;

    /**
     Integrates first-order GGA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed AO density matrix.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    auto _integrateVxcGradientForGGA(const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const std::vector<const double*>& rwDensityPointers,
                                     const std::vector<const double*>& gsDensityPointers,
                                     const CMolecularGrid&   molecularGrid,
                                     const CXCFunctional&    xcFunctional) const -> CDenseMatrix;

    auto _integrateVxcGradientForGGAOpenShell(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const std::vector<const double*>& rwDensityPointers,
                                              const std::vector<const double*>& gsDensityPointers,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCFunctional&    xcFunctional) const -> CDenseMatrix;

   public:
    /**
     Creates an XC integrator object.
     */
    CXCMolecularGradient();

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed AO density matrix (to
            be contracted with GTO gradient).
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const std::vector<const double*>& rwDensityPointers,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates first-order exchange-correlation functional contribution to
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
    CDenseMatrix integrateVxcPDFTGradient(const CMolecule&                    molecule,
                                          const CMolecularBasis&              basis,
                                          const double*                       densityMatrixPointer,
                                          const CDenseMatrix&                 twoBodyDensityMatrix,
                                          const CDenseMatrix&                 activeMOs,
                                          const CMolecularGrid&               molecularGrid,
                                          const CXCPairDensityFunctional&     fvxc,
                                          const double                        rs_omega) const;

};

#endif /* XCMolecularGradient_hpp */
