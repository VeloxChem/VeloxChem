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

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <mpi.h>

#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"

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

    /**
     Integrates first-order LDA exchange-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @return the AO Kohn-Sham matrix.
     */
    auto _integrateVxcFockForLDA(const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& densityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const std::string&      flag) const -> CAOKohnShamMatrix;

    auto _integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho) const -> CDenseMatrix;

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
     @param densityMatrix the AO density matrix object.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    auto integrateVxcFock(const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& densityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel) const -> CAOKohnShamMatrix;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    auto computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) const
        -> CDenseMatrix;
};

#endif /* XCIntegrator_hpp */
