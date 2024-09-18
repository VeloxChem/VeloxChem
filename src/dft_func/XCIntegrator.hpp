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

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <mpi.h>

#include <array>
#include <string>

#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

/**
 Class CXCIntegrator implements XC integrator.

 @author X. Li, K. Ahmadzadeh, M. Delcey
 */
class CXCIntegrator
{
    /**
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

   public:
    /**
     Creates an XC integrator object.

     @param comm the MPI communicator.
     */
    CXCIntegrator(MPI_Comm comm);

    /**
     Integrates first-order exchange-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to AO density matrices.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    auto integrateVxcFock(const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const std::string&                xcFuncLabel) const -> CAOKohnShamMatrix;

    /**
     Integrates first-order exchange-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to AO density matrices.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    auto integrateVxcFock(const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const CXCFunctional&              fvxc) const -> CAOKohnShamMatrix;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
     Fock matrix.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed density matrices.
     @param gsDensityPointers the pointers to ground state density matrices.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    auto integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                          const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const std::string&                xcFuncLabel) const -> void;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
     Fock matrix.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed density matrices.
     @param gsDensityPointers the pointers to ground state density matrices.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
    auto integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                          const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const CXCFunctional&              fvxc) const -> void;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    auto computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                      const CMolecularBasis& basis,
                                      const CMolecularGrid&  molecularGrid) const -> CDenseMatrix;
};

#endif /* XCIntegrator_hpp */