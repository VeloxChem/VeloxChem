//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularGrid.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

CMolecularGrid::CMolecularGrid()

    : _isDistributed(false)
{
}

CMolecularGrid::CMolecularGrid(const CMemBlock2D<double>& gridPoints)

    : _isDistributed(false)

    , _gridPoints(gridPoints)
{
}

CMolecularGrid::CMolecularGrid(const CMolecularGrid& source)

    : _isDistributed(source._isDistributed)

    , _gridPoints(source._gridPoints)
{
}

CMolecularGrid::CMolecularGrid(CMolecularGrid&& source) noexcept

    : _isDistributed(std::move(source._isDistributed))

    , _gridPoints(std::move(source._gridPoints))
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

    return *this;
}

CMolecularGrid&
CMolecularGrid::operator=(CMolecularGrid&& source) noexcept
{
    if (this == &source) return *this;

    _isDistributed = std::move(source._isDistributed);

    _gridPoints = std::move(source._gridPoints);

    return *this;
}

bool
CMolecularGrid::operator==(const CMolecularGrid& other) const
{
    if (_isDistributed != other._isDistributed) return false;

    if (_gridPoints != other._gridPoints) return false;

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
    if (nGridPoints < getNumberOfGridPoints())
    {
        _gridPoints = _gridPoints.slice(0, nGridPoints); 
    }
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
    if (!_isDistributed)
    {
        _isDistributed = true;

        _gridPoints.scatter(rank, nodes, comm);
    }
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
            std::fprintf(f, "%8i %16.12lf %16.12lf %16.12lf %16.12lf\n", i, rx[i], ry[i], rz[i], rw[i]);
        }
    }
    
    std::fclose(f);
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
