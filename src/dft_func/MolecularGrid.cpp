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

#include "MolecularGrid.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "GridPartitioner.hpp"

CMolecularGrid::CMolecularGrid()

    : _isDistributed(false)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(1024)
{
}

CMolecularGrid::CMolecularGrid(const int maxNumberOfGridPointsPerBox)

    : _isDistributed(false)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(maxNumberOfGridPointsPerBox)
{
}

CMolecularGrid::CMolecularGrid(const CDenseMatrix& gridPoints)

    : _isDistributed(false)

    , _gridPoints(gridPoints)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(1024)
{
}

CMolecularGrid::CMolecularGrid(const CDenseMatrix& gridPoints, const int maxNumberOfGridPointsPerBox)

    : _isDistributed(false)

    , _gridPoints(gridPoints)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(maxNumberOfGridPointsPerBox)
{
}

CMolecularGrid::CMolecularGrid(const CMolecularGrid& source)

    : _isDistributed(source._isDistributed)

    , _gridPoints(source._gridPoints)

    , _isPartitioned(source._isPartitioned)

    , _gridPointCounts(source._gridPointCounts)

    , _gridPointCountsOriginal(source._gridPointCountsOriginal)

    , _gridPointDisplacements(source._gridPointDisplacements)

    , _gridPointDisplacementsOriginal(source._gridPointDisplacementsOriginal)

    , _maxNumberOfGridPointsPerBox(source._maxNumberOfGridPointsPerBox)
{
}

CMolecularGrid::CMolecularGrid(CMolecularGrid&& source) noexcept

    : _isDistributed(std::move(source._isDistributed))

    , _gridPoints(std::move(source._gridPoints))

    , _isPartitioned(std::move(source._isPartitioned))

    , _gridPointCounts(std::move(source._gridPointCounts))

    , _gridPointCountsOriginal(std::move(source._gridPointCountsOriginal))

    , _gridPointDisplacements(std::move(source._gridPointDisplacements))

    , _gridPointDisplacementsOriginal(std::move(source._gridPointDisplacementsOriginal))

    , _maxNumberOfGridPointsPerBox(std::move(source._maxNumberOfGridPointsPerBox))
{
}

CMolecularGrid::~CMolecularGrid()
{
}

auto
CMolecularGrid::operator=(const CMolecularGrid& source) -> CMolecularGrid&
{
    if (this == &source) return *this;

    _isDistributed = source._isDistributed;

    _gridPoints = source._gridPoints;

    _isPartitioned = source._isPartitioned;

    _gridPointCounts = source._gridPointCounts;

    _gridPointCountsOriginal = source._gridPointCountsOriginal;

    _gridPointDisplacements = source._gridPointDisplacements;

    _gridPointDisplacementsOriginal = source._gridPointDisplacementsOriginal;

    _maxNumberOfGridPointsPerBox = source._maxNumberOfGridPointsPerBox;

    return *this;
}

auto
CMolecularGrid::operator=(CMolecularGrid&& source) noexcept -> CMolecularGrid&
{
    if (this == &source) return *this;

    _isDistributed = std::move(source._isDistributed);

    _gridPoints = std::move(source._gridPoints);

    _isPartitioned = std::move(source._isPartitioned);

    _gridPointCounts = std::move(source._gridPointCounts);

    _gridPointCountsOriginal = std::move(source._gridPointCountsOriginal);

    _gridPointDisplacements = std::move(source._gridPointDisplacements);

    _gridPointDisplacementsOriginal = std::move(source._gridPointDisplacementsOriginal);

    _maxNumberOfGridPointsPerBox = std::move(source._maxNumberOfGridPointsPerBox);

    return *this;
}

auto
CMolecularGrid::operator==(const CMolecularGrid& other) const -> bool
{
    if (_isDistributed != other._isDistributed) return false;

    if (_gridPoints != other._gridPoints) return false;

    if (_isPartitioned != other._isPartitioned) return false;

    if (_gridPointCounts != other._gridPointCounts) return false;

    if (_gridPointCountsOriginal != other._gridPointCountsOriginal) return false;

    if (_gridPointDisplacements != other._gridPointDisplacements) return false;

    if (_gridPointDisplacementsOriginal != other._gridPointDisplacementsOriginal) return false;

    if (_maxNumberOfGridPointsPerBox != other._maxNumberOfGridPointsPerBox) return false;

    return true;
}

auto
CMolecularGrid::operator!=(const CMolecularGrid& other) const -> bool
{
    return !(*this == other);
}

auto
CMolecularGrid::slice(const int nGridPoints) -> void
{
    std::string errpartitioned("MolecularGrid.slice: Cannot slice partitioned molecular grid");

    errors::assertMsgCritical(!_isPartitioned, errpartitioned);

    if (nGridPoints < getNumberOfGridPoints())
    {
        _gridPoints = _gridPoints.slice(0, nGridPoints);
    }
}

auto
CMolecularGrid::getGridPoints() const -> CDenseMatrix
{
    return _gridPoints;
}

auto
CMolecularGrid::getNumberOfGridPoints() const -> int
{
    return _gridPoints.getNumberOfColumns();
}

auto
CMolecularGrid::getCoordinatesX() const -> const double*
{
    return _gridPoints.row(0);
}

auto
CMolecularGrid::getCoordinatesY() const -> const double*
{
    return _gridPoints.row(1);
}

auto
CMolecularGrid::getCoordinatesZ() const -> const double*
{
    return _gridPoints.row(2);
}

auto
CMolecularGrid::getWeights() const -> const double*
{
    return _gridPoints.row(3);
}

auto
CMolecularGrid::getSpatialExtent() const -> std::array<double, 6>
{
    // initialize min, max coordinates

    double xmin = 0.0, ymin = 0.0, zmin = 0.0;

    double xmax = 0.0, ymax = 0.0, zmax = 0.0;

    auto ngpoints = getNumberOfGridPoints();

    if (ngpoints > 0)
    {
        // set up pointers to grid data

        auto rx = getCoordinatesX();

        auto ry = getCoordinatesY();

        auto rz = getCoordinatesZ();

        // update min, max coordinates

        xmin = rx[0];
        xmax = rx[0];

        ymin = ry[0];
        ymax = ry[0];

        zmin = rz[0];
        zmax = rz[0];

        for (int i = 0; i < ngpoints; i++)
        {
            if (xmin > rx[i]) xmin = rx[i];

            if (ymin > ry[i]) ymin = ry[i];

            if (zmin > rz[i]) zmin = rz[i];

            if (xmax < rx[i]) xmax = rx[i];

            if (ymax < ry[i]) ymax = ry[i];

            if (zmax < rz[i]) zmax = rz[i];
        }
    }

    std::array<double, 6> rext = {xmin, ymin, zmin, xmax, ymax, zmax};

    return rext;
}

auto
CMolecularGrid::partitionGridPoints() -> std::string
{
    if (!_isPartitioned)
    {
        _isPartitioned = true;

        auto boxdim = getSpatialExtent();

        // padding xmin, ymin & zmin

        for (int i = 0; i < 3; i++)
            boxdim[i] -= 0.0001;

        // padding xmax, ymax & zmax

        for (int i = 3; i < 6; i++)
            boxdim[i] += 0.0001;

        CGridPartitioner partitioner(CGridBox(boxdim, _gridPoints), _maxNumberOfGridPointsPerBox);

        partitioner.partitionGridPoints();

        _gridPoints = partitioner.getAllGridPoints();

        _gridPointCounts = partitioner.getGridPointCounts();

        _gridPointCountsOriginal = partitioner.getGridPointCounts();

        _gridPointDisplacements = partitioner.getGridPointDisplacements();

        _gridPointDisplacementsOriginal = partitioner.getGridPointDisplacements();

        return partitioner.getGridStatistics();
    }

    return std::string("");
}

auto
CMolecularGrid::distributeCountsAndDisplacements(const int rank, const int nnodes) -> void
{
    std::string errnotpartitioned("MolecularGrid.distributeCountsAndDisplacements: Cannot distribute unpartitioned molecular grid");

    errors::assertMsgCritical(_isPartitioned, errnotpartitioned);

    if (!_isDistributed)
    {
        _isDistributed = true;

        auto numboxes = static_cast<int>(_gridPointCountsOriginal.size());

        // sort before distribute

        std::vector<std::pair<int, int>> count_index_pairs;

        for (int box_id = 0; box_id < numboxes; box_id++)
        {
            auto count = _gridPointCountsOriginal[box_id];

            count_index_pairs.push_back(std::make_pair(count, box_id));
        }

        std::sort(count_index_pairs.begin(), count_index_pairs.end());

        // update original counts and displacements

        std::vector<int> orig_counts, orig_displs;

        for (int box_id = 0; box_id < numboxes; box_id++)
        {
            auto index = count_index_pairs[box_id].second;

            orig_counts.push_back(_gridPointCountsOriginal[index]);
            orig_displs.push_back(_gridPointDisplacementsOriginal[index]);
        }

        _gridPointCountsOriginal = orig_counts;
        _gridPointDisplacementsOriginal = orig_displs;

        // update counts and displacements

        std::vector<int> counts_p, displs_p;

        for (int box_id = numboxes - 1 - rank; box_id >= 0; box_id -= nnodes)
        {
            counts_p.push_back(_gridPointCountsOriginal[box_id]);
            displs_p.push_back(_gridPointDisplacementsOriginal[box_id]);
        }

        _gridPointCounts = counts_p;
        _gridPointDisplacements = displs_p;
    }
}

auto
CMolecularGrid::reDistributeCountsAndDisplacements(const int rank, const int nnodes) -> void
{
    std::string errpartitioned("MolecularGrid.reDistributeCountsAndDisplacements: Molecular grid must be parttioned and distributed");

    errors::assertMsgCritical((_isPartitioned && _isDistributed), errpartitioned);

    // recalculate counts and displacements

    auto numboxes = static_cast<int>(_gridPointCountsOriginal.size());

    std::vector<int> counts_p, displs_p;

    for (int box_id = numboxes - 1 - rank; box_id >= 0; box_id -= nnodes)
    {
        counts_p.push_back(_gridPointCountsOriginal[box_id]);
        displs_p.push_back(_gridPointDisplacementsOriginal[box_id]);
    }

    _gridPointCounts = counts_p;
    _gridPointDisplacements = displs_p;
}

auto
CMolecularGrid::isPartitioned() const -> bool
{
    return _isPartitioned;
}

auto
CMolecularGrid::getGridPointCounts() const -> std::vector<int>
{
    return _gridPointCounts;
}

auto
CMolecularGrid::getGridPointDisplacements() const -> std::vector<int>
{
    return _gridPointDisplacements;
}

auto
CMolecularGrid::getMaxNumberOfGridPointsPerBox() const -> int
{
    return _maxNumberOfGridPointsPerBox;
}
