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

#ifndef MolecularGrid_hpp
#define MolecularGrid_hpp

#include <array>
#include <cstdint>
#include <ostream>
#include <vector>

#include "DenseMatrix.hpp"

/**
 Class CMolecularGrid class generates molecular grid.

 @author Z. Rinkevicius
 */
class CMolecularGrid
{
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

    /**
     The displacement of grid points in the partitioned grid boxes.
     */
    std::vector<int64_t> _gridPointDisplacements;

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
     Creates a molecular grid object.

     @param gridPoints the 2D memory block object with grid points data.
     */
    CMolecularGrid(const CDenseMatrix& gridPoints);

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
    CMolecularGrid& operator=(const CMolecularGrid& source);

    /**
     Assigns a molecular grid object by moving other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid& operator=(CMolecularGrid&& source) noexcept;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are equal, false otherwise.
     */
    bool operator==(const CMolecularGrid& other) const;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularGrid& other) const;

    /**
     Reduces size of molecular grid by slicing all grid points beoynd given number of grid points.

     @param nGridPoints the number of grid points.
     */
    void slice(const int64_t nGridPoints);

    /**
     Gets grid points in molecular grid object.

     @return the grid points.
     */
    CDenseMatrix getGridPoints() const;

    /**
     Gets number of grid points in molecular grid object.

     @return the number of grid points.
     */
    int64_t getNumberOfGridPoints() const;

    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian X coordinates of grid points.
     */
    const double* getCoordinatesX() const;

    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.

     @return the pointer to Cartesian X coordinates of grid points.
     */
    double* getCoordinatesX();

    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.

     @return the constant  pointer to Cartesian Y coordinates of grid points.
     */
    const double* getCoordinatesY() const;

    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.

     @return the pointer to Cartesian Y coordinates of grid points.
     */
    double* getCoordinatesY();

    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian Z coordinates of grid points.
     */
    const double* getCoordinatesZ() const;

    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.

     @return the pointer to Cartesian Z coordinates of grid points.
     */
    double* getCoordinatesZ();

    /**
     Gets weights of grid points in molecular grid object.

     @return the constant pointer to weights of grid points.
     */
    const double* getWeights() const;

    /**
     Gets weights of grid points in molecular grid object.

     @return the pointer to weights of grid points.
     */
    double* getWeights();

    /**
     Gets spatial extent of molecular grid.

     @return the spatial extent (min x, min y, min z, max x, max y, max z).
     */
    std::array<double, 6> getSpatialExtent() const;

    /**
     Partitions grid points into boxes.

     @return summary of grid points partitioning as a string.
     */
    std::string partitionGridPoints();

    /**
     Checks whether the molecular grid has been partitioned.

     @return whether the molecular grid has been partitioned.
     */
    bool isPartitioned() const;

    /**
     Gets the grid point counts.

     @return the grid point counts.
     */
    std::vector<int64_t> getGridPointCounts() const;

    /**
     Gets the grid point displacements.

     @return the grid point displacements.
     */
    std::vector<int64_t> getGridPointDisplacements() const;

    /**
     Gets the maximum number of grid points in a grid box.

     @return the maximum number of grid points in a grid box.
     */
    int64_t getMaxNumberOfGridPointsPerBox() const;

    /**
     Converts molecular grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the molecular grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CMolecularGrid& source);
};

#endif /* MolecularGrid_hpp */
