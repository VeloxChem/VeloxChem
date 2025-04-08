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

#ifndef GridPartitioner_hpp
#define GridPartitioner_hpp

#include <list>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"
#include "GridBox.hpp"

/**
 Class CGridPartitioner implemention DFT grid points partitioner.
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
    int _numberOfPointsThreshold;

   public:
    /**
     Creates a grid partitioner object.
     */
    CGridPartitioner(const CGridBox& box, const int numGridPointsThreshold);

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
    auto findCenter(const double* ptr, const int num) const -> double;

    /**
     Finds the position for dividing the box.

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the one-dimensional median.
     */
    auto findDividingPosition(const double* ptr, const int num) const -> double;

    /**
     Estimates one-dimensional median using the P2 algorithm.
     https://doi.org/10.1145/4372.4378

     @param ptr the pointer to data.
     @param num the number of data points.
     @return the estimated one-dimensional median.
     */
    auto estimateMedianP2(const double* ptr, const int num) const -> double;

    /**
     Gets number of boxes.

     @return the number of boxes.
     */
    auto getNumberOfBoxes() const -> int;

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
    auto getTotalNumberOfGridPoints() const -> int;

    /**
     Gets all grid points in one object.

     @return all grid points.
     */
    auto getAllGridPoints() const -> CDenseMatrix;

    /**
     Gets grid point counts in the boxes.

     @return grid point counts.
     */
    auto getGridPointCounts() const -> std::vector<int>;

    /**
     Gets grid point displacements in the boxes.

     @return grid point displacements.
     */
    auto getGridPointDisplacements() const -> std::vector<int>;
};

#endif /* GridPartitioner_hpp */
