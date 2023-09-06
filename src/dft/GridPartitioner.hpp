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

#ifndef GridPartitioner_hpp
#define GridPartitioner_hpp

#include <list>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"
#include "GridBox.hpp"

/**
 Class CGridPartitioner implemention DFT grid points partitioner.

 @author X. Li
 */
class CGridPartitioner
{
    /**
     Grid boxes containing DFT grid points.
     */
    std::list<CGridBox> _boxes;

    /**
     Threshold for number of points in one box.
     */
    int64_t _numberOfPointsThreshold;

   public:
    /**
     Creates a grid partitioner object.
     */
    CGridPartitioner(const CGridBox& box, const int64_t numGridPointsThreshold);

    /**
     Partitions the grid points.
     */
    auto partitionGridPoints() -> void;

    /**
     Divides a grid box into two boxes.

     @param box the grid box.
     @return a list of smaller boxes.
     */
    auto bisectBox(const CGridBox& box) const -> std::list<CGridBox>;

    /**
     Finds one-dimensional center of points.

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the one-dimensional center.
     */
    auto findCenter(const double* ptr, const int64_t num) const -> double;

    /**
     Finds the position for dividing the box.

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the one-dimensional median.
     */
    auto findDividingPosition(const double* ptr, const int64_t num) const -> double;

    /**
     Estimates one-dimensional median using the P2 algorithm.
     https://doi.org/10.1145/4372.4378

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the estimated one-dimensional median.
     */
    auto estimateMedianP2(const double* ptr, const int64_t num) const -> double;

    /**
     Gets number of boxes.

     @return the number of boxes.
     */
    auto getNumberOfBoxes() const -> int64_t;

    /**
     Gets information about grid boxes.

     @return a string containing information about grid boxes.
     */
    auto getGridInformation() const -> std::string;

    /**
     Gets statistics about grid boxes.

     @return a string containing statistics about grid boxes.
     */
    auto getGridStatistics() const -> std::string;

    /**
     Gets the total number of grid points.

     @return the total number of grid points.
     */
    auto getTotalNumberOfGridPoints() const -> int64_t;

    /**
     Gets all grid points in one object.

     @return all grid points.
     */
    auto getAllGridPoints() const -> CDenseMatrix;

    /**
     Gets grid point counts in the boxes.

     @return grid point counts.
     */
    auto getGridPointCounts() const -> std::vector<int64_t>;

    /**
     Gets grid point displacements in the boxes.

     @return grid point displacements.
     */
    auto getGridPointDisplacements() const -> std::vector<int64_t>;
};

#endif /* GridPartitioner_hpp */
