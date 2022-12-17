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
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <mpi.h>

#include "DenseLinearAlgebra.hpp"
#include "ErrorHandler.hpp"
#include "GridPartitioner.hpp"
#include "GtoContainer.hpp"
#include "GtoFunc.hpp"
#include "XCVarsType.hpp"

CMolecularGrid::CMolecularGrid()

    : _isDistributed(false)

    , _isPartitioned(false)

    , _maxNumberOfGridPointsPerBox(1024)
{
}

CMolecularGrid::CMolecularGrid(const CMemBlock2D<double>& gridPoints)

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
CMolecularGrid::slice(const int32_t nGridPoints)
{
    std::string errpartitioned("MolecularGrid.slice: Cannot slice partitioned molecular grid");

    errors::assertMsgCritical(!_isPartitioned, errpartitioned);

    if (nGridPoints < getNumberOfGridPoints())
    {
        _gridPoints = _gridPoints.slice(0, nGridPoints); 
    }
}

CMemBlock2D<double>
CMolecularGrid::getGridPoints() const
{
    return _gridPoints;
}

int32_t
CMolecularGrid::getNumberOfGridPoints() const
{
    return _gridPoints.size(0);
}

const double*
CMolecularGrid::getCoordinatesX() const
{
    return _gridPoints.data(0);
}

double*
CMolecularGrid::getCoordinatesX()
{
    return _gridPoints.data(0);
}

const double*
CMolecularGrid::getCoordinatesY() const
{
    return _gridPoints.data(1);
}

double*
CMolecularGrid::getCoordinatesY()
{
    return _gridPoints.data(1);
}

const double*
CMolecularGrid::getCoordinatesZ() const
{
    return _gridPoints.data(2);
}

double*
CMolecularGrid::getCoordinatesZ()
{
    return _gridPoints.data(2);
}

const double*
CMolecularGrid::getWeights() const
{
    return _gridPoints.data(3);
}

double*
CMolecularGrid::getWeights()
{
    return _gridPoints.data(3);
}

void
CMolecularGrid::distribute(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    std::string errpartitioned("MolecularGrid.distribute: Cannot distribute partitioned molecular grid");

    errors::assertMsgCritical(!_isPartitioned, errpartitioned);

    if (!_isDistributed)
    {
        _isDistributed = true;

        _gridPoints.scatter(rank, nodes, comm);
    }
}

void
CMolecularGrid::broadcast(int32_t rank, MPI_Comm comm)
{
    std::string errpartitioned("MolecularGrid.broadcast: Cannot broadcast partitioned molecular grid");

    errors::assertMsgCritical(!_isPartitioned, errpartitioned);

    _gridPoints.broadcast(rank, comm);
}

void
CMolecularGrid::read_blocked_grid(const std::string& fileName)
{
    const int32_t mpoints = 20000000;
    
    CMemBlock2D<double> rgrid(mpoints, 4);
    
    auto f = std::fopen(fileName.c_str(), "rb");
    
    if (f != nullptr)
    {
        int32_t npnt = 0;
        
        int32_t cpnt = 0;
        
        while(true)
        {
            // read number of grid points in block
            
            auto nblock = std::fread(&npnt, sizeof(int32_t), 1u, f);
            
            if (nblock < 1u) break;
            
            // read shell block data
            
            int32_t nshells = 0;
            
            nblock = std::fread(&nshells, sizeof(int32_t), 1u, f);
            
            if (nblock < 1u) break;
            
            for(int32_t i = 0; i < 2 * nshells; i++)
            {
                int32_t bf = 0;
                
                nblock = std::fread(&bf, sizeof(int32_t), 1u, f);
                
                if (nblock < 1u) break;
            }
            
            // read coordinates and weights

            auto rx = rgrid.data(0, cpnt);
            
            auto ry = rgrid.data(1, cpnt);
            
            auto rz = rgrid.data(2, cpnt);
            
            auto rw = rgrid.data(3, cpnt);
            
            nblock = std::fread(rx, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(ry, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(rz, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(rw, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            cpnt += npnt;
        }
        
        _isDistributed = false;
        
        _gridPoints = rgrid.slice(0, cpnt);
    }
    
    std::fclose(f);
}

void
CMolecularGrid::read_raw_grid(const std::string& fileName)
{
    // allocate raw grid
    
    const int32_t mpoints = 20000000;
    
    CMemBlock2D<double> rgrid(mpoints, 4);
    
    // set up pointers to grid data
    
    auto rx = rgrid.data(0);
    
    auto ry = rgrid.data(1);
    
    auto rz = rgrid.data(2);
    
    auto rw = rgrid.data(3);
    
    // read grid data from text file
    
    std::ifstream fst(fileName.c_str());
    
    int32_t npoints = 0;
    
    while(true)
    {
        std::string str;
        
        std::getline(fst, str);

        if (fst.eof()) break;
        
        std::istringstream iss(str);
        
        int32_t idx = 0;
        
        iss >> idx >> rx[npoints] >> ry[npoints] >> rz[npoints] >> rw[npoints];
        
        npoints++;
    }
    
    _isDistributed = false;
    
    _gridPoints = rgrid.slice(0, npoints);
}

void
CMolecularGrid::write_raw_grid(const std::string fileName) const
{
    auto f = std::fopen(fileName.c_str(), "w");
    
    if (f != nullptr)
    {
        auto rx = getCoordinatesX();
        
        auto ry = getCoordinatesY();
        
        auto rz = getCoordinatesZ();
        
        auto rw = getWeights();
        
        for (int32_t i = 0; i < getNumberOfGridPoints(); i++)
        {
            std::fprintf(f, "%8i %20.16lf %20.16lf %20.16lf %20.16lf\n", i, rx[i], ry[i], rz[i], rw[i]);
        }
    }
    
    std::fclose(f);
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
        
        xmin = rx[0]; xmax = rx[0];
        
        ymin = ry[0]; ymax = ry[0];
        
        zmin = rz[0]; zmax = rz[0];
        
        for (int32_t i = 0; i < ngpoints; i++)
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

        for (int32_t i = 0; i < 3; i++) boxdim[i] -= 0.0001;

        // padding xmax, ymax & zmax

        for (int32_t i = 3; i < 6; i++) boxdim[i] += 0.0001;

        CGridPartitioner partitioner(CGridBox(boxdim, _gridPoints), _maxNumberOfGridPointsPerBox);

        partitioner.partitionGridPoints();

        _gridPoints = partitioner.getAllGridPoints();

        _gridPointCounts = partitioner.getGridPointCounts();

        _gridPointDisplacements = partitioner.getGridPointDisplacements();

        return partitioner.getGridStatistics();
    }

    return std::string("");
}

void
CMolecularGrid::distributeCountsAndDisplacements(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    std::string errnotpartitioned("MolecularGrid.distributeCountsAndDisplacements: Cannot broadcast unpartitioned molecular grid");

    errors::assertMsgCritical(_isPartitioned, errnotpartitioned);

    if (!_isDistributed)
    {
        _isDistributed = true;

        auto numboxes = _gridPointCounts.size();

        // sort before distribute

        std::vector<std::pair<int32_t, int32_t>> count_index_pairs;

        for (int32_t box_id = 0; box_id < numboxes; box_id++)
        {
            auto count = _gridPointCounts.data()[box_id];

            count_index_pairs.push_back(std::make_pair(count, box_id));
        }

        std::sort(count_index_pairs.begin(), count_index_pairs.end());

        // re-arrange sorted indices for parallelization

        std::vector<int32_t> newcounts;

        std::vector<int32_t> newdispls;

        for (int32_t p = 0; p < nodes; p++)
        {
            for (int32_t box_id = numboxes - 1 - p; box_id >= 0; box_id -= nodes)
            {
                auto index = count_index_pairs[box_id].second;

                newcounts.push_back(_gridPointCounts.data()[index]);

                newdispls.push_back(_gridPointDisplacements.data()[index]);
            }
        }

        std::memcpy(_gridPointCounts.data(), newcounts.data(), numboxes * sizeof(int32_t));

        std::memcpy(_gridPointDisplacements.data(), newdispls.data(), numboxes * sizeof(int32_t));

        // distribute counts and displacements

        _gridPointCounts.scatter(rank, nodes, comm);

        _gridPointDisplacements.scatter(rank, nodes, comm);

        // broadcast all grid points, since counts and displacements are already distributed

        _gridPoints.broadcast(rank, comm);
    }
}

void
CMolecularGrid::reDistributeCountsAndDisplacements(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    _isPartitioned = false;

    partitionGridPoints();

    _isDistributed = false;

    distributeCountsAndDisplacements(rank, nodes, comm);
}

bool
CMolecularGrid::isPartitioned() const
{
    return _isPartitioned;
}

CMemBlock<int32_t>
CMolecularGrid::getGridPointCounts() const
{
    return _gridPointCounts;
}

CMemBlock<int32_t>
CMolecularGrid::getGridPointDisplacements() const
{
    return _gridPointDisplacements;
}

int32_t
CMolecularGrid::getMaxNumberOfGridPointsPerBox() const
{
    return _maxNumberOfGridPointsPerBox;
}

std::ostream&
operator<<(std::ostream& output, const CMolecularGrid& source)
{
    output << std::endl;
    
    output << "[CMolecularGrid (Object):" << &source << "]" << std::endl;

    output << "_isDistributed: " << source._isDistributed << std::endl;

    output << "_gridPoints: " << std::endl;

    output << source._gridPoints << std::endl;

    return output;
}
