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

#include "MolecularGrid.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
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

CMolecularGrid::CMolecularGrid(const CDenseMatrix& gridPoints)

    : _isDistributed(false)

    , _gridPoints(gridPoints)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(1024)
{
}

CMolecularGrid::CMolecularGrid(const CMolecularGrid& source)

    : _isDistributed(source._isDistributed)

    , _gridPoints(source._gridPoints)

    , _isPartitioned(source._isPartitioned)

    , _gridPointCounts(source._gridPointCounts)

    , _gridPointDisplacements(source._gridPointDisplacements)

    , _maxNumberOfGridPointsPerBox(source._maxNumberOfGridPointsPerBox)
{
}

CMolecularGrid::CMolecularGrid(CMolecularGrid&& source) noexcept

    : _isDistributed(std::move(source._isDistributed))

    , _gridPoints(std::move(source._gridPoints))

    , _isPartitioned(std::move(source._isPartitioned))

    , _gridPointCounts(std::move(source._gridPointCounts))

    , _gridPointDisplacements(std::move(source._gridPointDisplacements))

    , _maxNumberOfGridPointsPerBox(std::move(source._maxNumberOfGridPointsPerBox))
{
}

CMolecularGrid::~CMolecularGrid()
{
}

CMolecularGrid&
CMolecularGrid::operator=(const CMolecularGrid& source)
{
    if (this == &source) return *this;

    _isDistributed = source._isDistributed;

    _gridPoints = source._gridPoints;

    _isPartitioned = source._isPartitioned;

    _gridPointCounts = source._gridPointCounts;

    _gridPointDisplacements = source._gridPointDisplacements;

    _maxNumberOfGridPointsPerBox = source._maxNumberOfGridPointsPerBox;

    return *this;
}

CMolecularGrid&
CMolecularGrid::operator=(CMolecularGrid&& source) noexcept
{
    if (this == &source) return *this;

    _isDistributed = std::move(source._isDistributed);

    _gridPoints = std::move(source._gridPoints);

    _isPartitioned = std::move(source._isPartitioned);

    _gridPointCounts = std::move(source._gridPointCounts);

    _gridPointDisplacements = std::move(source._gridPointDisplacements);

    _maxNumberOfGridPointsPerBox = std::move(source._maxNumberOfGridPointsPerBox);

    return *this;
}

bool
CMolecularGrid::operator==(const CMolecularGrid& other) const
{
    if (_isDistributed != other._isDistributed) return false;

    if (_gridPoints != other._gridPoints) return false;

    if (_isPartitioned != other._isPartitioned) return false;

    if (_gridPointCounts != other._gridPointCounts) return false;

    if (_gridPointDisplacements != other._gridPointDisplacements) return false;

    if (_maxNumberOfGridPointsPerBox != other._maxNumberOfGridPointsPerBox) return false;

    return true;
}

bool
CMolecularGrid::operator!=(const CMolecularGrid& other) const
{
    return !(*this == other);
}

void
CMolecularGrid::slice(const int64_t nGridPoints)
{
    std::string errpartitioned("MolecularGrid.slice: Cannot slice partitioned molecular grid");

    errors::assertMsgCritical(!_isPartitioned, errpartitioned);

    if (nGridPoints < getNumberOfGridPoints())
    {
        _gridPoints = _gridPoints.slice(0, nGridPoints);
    }
}

CDenseMatrix
CMolecularGrid::getGridPoints() const
{
    return _gridPoints;
}

int64_t
CMolecularGrid::getNumberOfGridPoints() const
{
    return _gridPoints.getNumberOfColumns();
}

const double*
CMolecularGrid::getCoordinatesX() const
{
    return _gridPoints.row(0);
}

double*
CMolecularGrid::getCoordinatesX()
{
    return _gridPoints.row(0);
}

const double*
CMolecularGrid::getCoordinatesY() const
{
    return _gridPoints.row(1);
}

double*
CMolecularGrid::getCoordinatesY()
{
    return _gridPoints.row(1);
}

const double*
CMolecularGrid::getCoordinatesZ() const
{
    return _gridPoints.row(2);
}

double*
CMolecularGrid::getCoordinatesZ()
{
    return _gridPoints.row(2);
}

const double*
CMolecularGrid::getWeights() const
{
    return _gridPoints.row(3);
}

double*
CMolecularGrid::getWeights()
{
    return _gridPoints.row(3);
}

std::array<double, 6>
CMolecularGrid::getSpatialExtent() const
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

        for (int64_t i = 0; i < ngpoints; i++)
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

std::string
CMolecularGrid::partitionGridPoints()
{
    if (!_isPartitioned)
    {
        _isPartitioned = true;

        auto boxdim = getSpatialExtent();

        // padding xmin, ymin & zmin

        for (int64_t i = 0; i < 3; i++)
            boxdim[i] -= 0.0001;

        // padding xmax, ymax & zmax

        for (int64_t i = 3; i < 6; i++)
            boxdim[i] += 0.0001;

        CGridPartitioner partitioner(CGridBox(boxdim, _gridPoints), _maxNumberOfGridPointsPerBox);

        partitioner.partitionGridPoints();

        _gridPoints = partitioner.getAllGridPoints();

        _gridPointCounts = partitioner.getGridPointCounts();

        _gridPointDisplacements = partitioner.getGridPointDisplacements();

        return partitioner.getGridStatistics();
    }

    return std::string("");
}

bool
CMolecularGrid::isPartitioned() const
{
    return _isPartitioned;
}

std::vector<int64_t>
CMolecularGrid::getGridPointCounts() const
{
    return _gridPointCounts;
}

std::vector<int64_t>
CMolecularGrid::getGridPointDisplacements() const
{
    return _gridPointDisplacements;
}

int64_t
CMolecularGrid::getMaxNumberOfGridPointsPerBox() const
{
    return _maxNumberOfGridPointsPerBox;
}

std::ostream&
operator<<(std::ostream& output, const CMolecularGrid& source)
{
    output << std::endl;

    output << "[CMolecularGrid (Object):" << &source << "]" << std::endl;

    output << "_isPartitioned: " << source._isPartitioned << std::endl;

    output << "_isDistributed: " << source._isDistributed << std::endl;

    output << "_gridPoints: " << std::endl;

    output << source._gridPoints << std::endl;

    return output;
}
