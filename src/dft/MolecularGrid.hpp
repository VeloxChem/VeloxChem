//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    std::vector<int64_t> _gridPointCounts;

    std::vector<int64_t> _gridPointCountsOriginal;

    /**
     The displacement of grid points in the partitioned grid boxes.
     */
    std::vector<int64_t> _gridPointDisplacements;

    std::vector<int64_t> _gridPointDisplacementsOriginal;

    /**
     The maximum number of grid points in a grid box.
     */
    int64_t _maxNumberOfGridPointsPerBox;

   public:
    /**
     Creates an empty molecular grid object.
     */
    CMolecularGrid();

    /**
     Creates an empty molecular grid object.

     @param maxNumberOfGridPointsPerBox the maximum number of grid points per box.
     */
    CMolecularGrid(const int64_t maxNumberOfGridPointsPerBox);

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
    CMolecularGrid(const CDenseMatrix& gridPoints, const int64_t maxNumberOfGridPointsPerBox);

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
    auto slice(const int64_t nGridPoints) -> void;

    /**
     Gets grid points in molecular grid object.

     @return the grid points.
     */
    auto getGridPoints() const -> CDenseMatrix;

    /**
     Gets number of grid points in molecular grid object.

     @return the number of grid points.
     */
    auto getNumberOfGridPoints() const -> int64_t;

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
    auto distributeCountsAndDisplacements(const int64_t rank, const int64_t nnodes) -> void;

    /**
     Redo distributing grid point counts and displacements within domain of MPI
     communacator and sets distribution flag to true.

     @param rank the MPI rank.
     @param nnodes the number of MPI processes.
     */
    auto reDistributeCountsAndDisplacements(const int64_t rank, const int64_t nnodes) -> void;

    /**
     Checks whether the molecular grid has been partitioned.

     @return whether the molecular grid has been partitioned.
     */
    auto isPartitioned() const -> bool;

    /**
     Gets the grid point counts.

     @return the grid point counts.
     */
    auto getGridPointCounts() const -> std::vector<int64_t>;

    /**
     Gets the grid point displacements.

     @return the grid point displacements.
     */
    auto getGridPointDisplacements() const -> std::vector<int64_t>;

    /**
     Gets the maximum number of grid points in a grid box.

     @return the maximum number of grid points in a grid box.
     */
    auto getMaxNumberOfGridPointsPerBox() const -> int64_t;
};

#endif /* MolecularGrid_hpp */
