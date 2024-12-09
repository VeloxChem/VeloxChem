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

#ifndef MolecularGrid_hpp
#define MolecularGrid_hpp

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"

/**
 Class CMolecularGrid class generates molecular grid.
 */
class CMolecularGrid
{
   private:
    /**
     The distribution status of molecular grid object.
     */
    bool _isDistributed;

    /**
     The grid points (coordinates, weights).
     */
    CDenseMatrix _gridPoints;

    /**
     The parition status of molecular grid object.
     */
    bool _isPartitioned;

    /**
     The number of grid points in the partitioned grid boxes.
     */
    std::vector<int> _gridPointCounts;

    /**
     The displacement of grid points in the partitioned grid boxes.
     */
    std::vector<int> _gridPointDisplacements;

    /**
     The maximum number of grid points in a grid box.
     */
    int _maxNumberOfGridPointsPerBox;

   public:
    /**
     Creates an empty molecular grid object.
     */
    CMolecularGrid();

    /**
     Creates an empty molecular grid object.

     @param maxNumberOfGridPointsPerBox the maximum number of grid points per box.
     */
    CMolecularGrid(const int maxNumberOfGridPointsPerBox);

    /**
     Creates a molecular grid object.

     @param gridPoints the 2D memory block object with grid points data.
     */
    CMolecularGrid(const CDenseMatrix& gridPoints);

    /**
     Creates a molecular grid object.

     @param gridPoints the 2D memory block object with grid points data.
     @param maxNumberOfGridPointsPerBox the maximum number of grid points per box.
     */
    CMolecularGrid(const CDenseMatrix& gridPoints, const int maxNumberOfGridPointsPerBox);

    /**
     Creates a molecular grid object by copying other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid(const CMolecularGrid& source);

    /**
     Creates a molecular grid object by by moving other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid(CMolecularGrid&& source) noexcept;

    /**
     Destroys a molecular grid object.
     */
    ~CMolecularGrid();

    /**
     Assigns a molecular grid object by copying other molecular grid object.

     @param source the molecular grid object.
     */
    auto operator=(const CMolecularGrid& source) -> CMolecularGrid&;

    /**
     Assigns a molecular grid object by moving other molecular grid object.

     @param source the molecular grid object.
     */
    auto operator=(CMolecularGrid&& source) noexcept -> CMolecularGrid&;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are equal, false otherwise.
     */
    auto operator==(const CMolecularGrid& other) const -> bool;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are not equal, false otherwise.
     */
    auto operator!=(const CMolecularGrid& other) const -> bool;

    /**
     Reduces size of molecular grid by slicing all grid points beoynd given number of grid points.

     @param nGridPoints the number of grid points.
     */
    auto slice(const int nGridPoints) -> void;

    /**
     Gets grid points in molecular grid object.

     @return the grid points.
     */
    auto getGridPoints() const -> CDenseMatrix;

    /**
     Gets number of grid points in molecular grid object.

     @return the number of grid points.
     */
    auto getNumberOfGridPoints() const -> int;

    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian X coordinates of grid points.
     */
    auto getCoordinatesX() const -> const double*;

    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.

     @return the constant  pointer to Cartesian Y coordinates of grid points.
     */
    auto getCoordinatesY() const -> const double*;

    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian Z coordinates of grid points.
     */
    auto getCoordinatesZ() const -> const double*;

    /**
     Gets weights of grid points in molecular grid object.

     @return the constant pointer to weights of grid points.
     */
    auto getWeights() const -> const double*;

    /**
     Gets spatial extent of molecular grid.

     @return the spatial extent (min x, min y, min z, max x, max y, max z).
     */
    auto getSpatialExtent() const -> std::array<double, 6>;

    /**
     Partitions grid points into boxes.

     @return summary of grid points partitioning as a string.
     */
    auto partitionGridPoints() -> std::string;

    /**
     Distributes grid point counts and displacements within domain of MPI
     communacator and sets distribution flag to true.

     @param rank the MPI rank.
     @param nnodes the number of MPI processes.
     */
    auto distributeCountsAndDisplacements(const int rank, const int nnodes) -> void;

    /**
     Redo distributing grid point counts and displacements within domain of MPI
     communacator and sets distribution flag to true.

     @param rank the MPI rank.
     @param nnodes the number of MPI processes.
     */
    auto reDistributeCountsAndDisplacements(const int rank, const int nnodes) -> void;

    /**
     Checks whether the molecular grid has been partitioned.

     @return whether the molecular grid has been partitioned.
     */
    auto isPartitioned() const -> bool;

    /**
     Gets the grid point counts.

     @return the grid point counts.
     */
    auto getGridPointCounts() const -> std::vector<int>;

    /**
     Gets the grid point displacements.

     @return the grid point displacements.
     */
    auto getGridPointDisplacements() const -> std::vector<int>;

    /**
     Gets the maximum number of grid points in a grid box.

     @return the maximum number of grid points in a grid box.
     */
    auto getMaxNumberOfGridPointsPerBox() const -> int;
};

#endif /* MolecularGrid_hpp */
