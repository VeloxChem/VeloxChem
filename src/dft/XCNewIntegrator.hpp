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

#ifndef XCNewIntegrator_hpp
#define XCNewIntegrator_hpp

#include <list>
#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "XCGradientGrid.hpp"

/**
 Class CXCNewIntegrator implements XC integrator.

 @author X. Li
 */
class CXCNewIntegrator
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     Grid boxes containing DFT grid points.
     */
    std::list<CGridBox> _boxes;

    /**
     Threshold for number of points in one box.
     */
    int32_t _numberOfPointsThreshold;

    /**
     Generates density grid.

     @param npoints the number of grid points.
     @param gtoValuesOnGridPoints the GTO values on grid points.
     @param densityMatrix the AO density matrix.
     @param xcFunType the type of exchange-correlation functional.
     @return the density grid.
     */
    CDensityGrid _generateDensityGrid(const int32_t           npoints,
                                      const CDenseMatrix&     gtoValuesOnGridPoints,
                                      const CAODensityMatrix& densityMatrix,
                                      const xcfun             xcFunType) const;

    /**
     Integrates first-order exchnage-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param xcoords the X coordinates of grid points.
     @param ycoords the Y coordinates of grid points.
     @param zcoords the Z coordinates of grid points.
     @param weights the weights of grid points.
     @param gtoValuesOnGridPoints the label of exchange-correlation functional.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForLDA(const int32_t          npoints,
                                                const double*          xcoords,
                                                const double*          ycoords,
                                                const double*          zcoords,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValuesOnGridPoints,
                                                const CXCGradientGrid& xcGradientGrid) const;

   public:
    /**
     Creates an XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCNewIntegrator(MPI_Comm comm);

    /**
     Destroys an XC integrator object.
     */
    ~CXCNewIntegrator();

    /**
     Sets threshold for number of points in one box.

     @param thresh the threshold for number of points in one box
     */
    void setNumberOfPointsThreshold(const int32_t thresh);

    /**
     Gets number of boxes.

     @return the number of boxes.
     */
    int32_t getNumberOfBoxes() const;

    /**
     Initializes molecular grid.

     @param molgrid the molecular grid.
     */
    void initializeGrid(const CMolecularGrid& molgrid);

    /**
     Partitions the grid points.
     */
    void partitionGrid();

    /**
     Divides a grid box into eight boxes.

     @param box the grid box.
     @return a list of smaller boxes.
     */
    std::list<CGridBox> divideBoxIntoEight(const CGridBox& box) const;

    /**
     Divides a grid box into two boxes.

     @param box the grid box.
     @return a list of smaller boxes.
     */
    std::list<CGridBox> divideBoxIntoTwo(const CGridBox& box) const;

    /**
     Gets information about grid boxes.

     @return a string containing information about grid boxes.
     */
    std::string getGridInformation() const;

    /**
     Gets statistics about grid boxes.

     @return a string containing statistics about grid boxes.
     */
    std::string getGridStatistics() const;

    /**
     Integrates first-order exchnage-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix integrateVxcFock(const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const CAODensityMatrix& densityMatrix,
                                       const std::string&      xcFuncLabel) const;
};

#endif /* XCNewIntegrator_hpp */
