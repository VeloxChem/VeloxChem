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

#include <vector>
#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

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
     Boxes (xmin,ymin,zmin,xmax,ymax,zmax) containing grid points.
     */
    std::vector<std::array<double, 6>> _boxes;

    /**
     Grid points (x,y,z,w) in the boxes.
     */
    std::vector<CMemBlock2D<double>> _pointsInBoxes;

    /**
     Threshold for number of points in one box.
     */
    int32_t _numberOfPointsThreshold;

    /**
     Maximum level for dividing the boxes.
     */
    int32_t _maximumLevel;

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
     Sets the maximum level of dividing the boxes.

     @param maxlevel the maximum level.
     */
    void setMaximumLevel(const int32_t maxlevel);

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

    void partitionGrid();

    std::string getGridInformation() const;

    CDenseMatrix integrateVxc(const CMolecule&        molecule,
                              const CMolecularBasis&  basis,
                              const CAODensityMatrix& densityMatrix,
                              const std::string&      xcFuncLabel) const;
};

#endif /* XCNewIntegrator_hpp */
