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

#include "GridBox.hpp"
#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"

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
    int32_t _numberOfPointsThreshold;

   public:
    /**
     Creates a grid partitioner object.
     */
    CGridPartitioner(const CGridBox& box);

    /**
     Destroys a grid partitioner object.
     */
    ~CGridPartitioner();

    /**
     Partitions the grid points.
     */
    void partitionGridPoints();

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
     Finds one-dimensional center of points.

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the one-dimensional center.
     */
    double findCenter(const double* ptr, const int32_t num) const;

    /**
     Finds one-dimensional median.

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the one-dimensional median.
     */
    double findMedian(const double* ptr, const int32_t num) const;

    /**
     Estimates one-dimensional median using the P2 algorithm.
     https://doi.org/10.1145/4372.4378

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the estimated one-dimensional median.
     */
    double estimateMedianP2(const double* ptr, const int32_t num) const;

    /**
     Gets number of boxes.

     @return the number of boxes.
     */
    int32_t getNumberOfBoxes() const;

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
     Gets the total number of grid points.

     @return the total number of grid points.
     */
    int32_t getTotalNumberOfGridPoints() const;

    /**
     Gets all grid points in one object.

     @return all grid points.
     */
    CMemBlock2D<double> getAllGridPoints() const;

    /**
     Gets grid point counts in the boxes.

     @return grid point counts.
     */
    CMemBlock<int32_t> getGridPointCounts() const;

    /**
     Gets grid point displacements in the boxes.

     @return grid point displacements.
     */
    CMemBlock<int32_t> getGridPointDisplacements() const;
};

#endif /* GridPartitioner_hpp */
