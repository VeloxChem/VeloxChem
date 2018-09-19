//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularGrid.hpp"

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

const double*
CMolecularGrid::getCoordinatesY() const
{
    return _gridPoints.data(1);
}

const double*
CMolecularGrid::getCoordinatesZ() const
{
    return _gridPoints.data(2);
}

const double*
CMolecularGrid::getWeights() const
{
    return _gridPoints.data(3);
}

void
CMolecularGrid::distribute(int32_t  rank,
                           int32_t  nodes,
                           MPI_Comm comm)
{
    if (!_isDistributed)
    {
        _isDistributed = true;

        _gridPoints.scatter(rank, nodes, comm);
    }
}

std::ostream&
operator<<(      std::ostream&   output,
           const CMolecularGrid& source)
{
    output << std::endl;

    output << "_isDistributed: " << source._isDistributed << std::endl;

    output << "_gridPoints: " << std::endl;

    output << source._gridPoints << std::endl;

    return output;
}
