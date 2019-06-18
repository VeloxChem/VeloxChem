//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CubicGrid.hpp"

#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MemBlock.hpp"

CCubicGrid::CCubicGrid()

    : _origin(std::vector<double>({0.0, 0.0, 0.0}))

    , _stepSize(std::vector<double>({0.0, 0.0, 0.0}))

    , _numPoints(std::vector<int32_t>({0, 0, 0}))

    , _values(CMemBlock<double>())
{
}

CCubicGrid::CCubicGrid(const std::vector<double>& origin, const std::vector<double>& stepSize, const std::vector<int32_t>& numPoints)

    : _origin(origin)

    , _stepSize(stepSize)

    , _numPoints(numPoints)
{
    std::string errmsg("CubicGrid: incorrect dimension");

    errors::assertMsgCritical(_origin.size() == 3, errmsg);

    errors::assertMsgCritical(_stepSize.size() == 3, errmsg);

    errors::assertMsgCritical(_numPoints.size() == 3, errmsg);

    _values = CMemBlock<double>(numPoints[0] * numPoints[1] * numPoints[2]);
}

CCubicGrid::CCubicGrid(const CCubicGrid& source)

    : _origin(source._origin)

    , _stepSize(source._stepSize)

    , _numPoints(source._numPoints)

    , _values(source._values)
{
}

CCubicGrid::CCubicGrid(CCubicGrid&& source) noexcept

    : _origin(std::move(source._origin))

    , _stepSize(std::move(source._stepSize))

    , _numPoints(std::move(source._numPoints))

    , _values(std::move(source._values))
{
}

CCubicGrid::~CCubicGrid()
{
}

CCubicGrid&
CCubicGrid::operator=(const CCubicGrid& source)
{
    if (this == &source) return *this;

    _origin = source._origin;

    _stepSize = source._stepSize;

    _numPoints = source._numPoints;

    _values = source._values;

    return *this;
}

CCubicGrid&
CCubicGrid::operator=(CCubicGrid&& source) noexcept
{
    if (this == &source) return *this;

    _origin = std::move(source._origin);

    _stepSize = std::move(source._stepSize);

    _numPoints = std::move(source._numPoints);

    _values = std::move(source._values);

    return *this;
}

bool
CCubicGrid::operator==(const CCubicGrid& other) const
{
    if (_origin != other._origin) return false;

    if (_stepSize != other._stepSize) return false;

    if (_numPoints != other._numPoints) return false;

    if (_values != other._values) return false;

    return true;
}

bool
CCubicGrid::operator!=(const CCubicGrid& other) const
{
    return !(*this == other);
}

double
CCubicGrid::originX() const
{
    return _origin[0];
}

double
CCubicGrid::originY() const
{
    return _origin[1];
}

double
CCubicGrid::originZ() const
{
    return _origin[2];
}

double
CCubicGrid::stepSizeX() const
{
    return _stepSize[0];
}

double
CCubicGrid::stepSizeY() const
{
    return _stepSize[1];
}

double
CCubicGrid::stepSizeZ() const
{
    return _stepSize[2];
}

int32_t
CCubicGrid::numPointsX() const
{
    return _numPoints[0];
}

int32_t
CCubicGrid::numPointsY() const
{
    return _numPoints[1];
}

int32_t
CCubicGrid::numPointsZ() const
{
    return _numPoints[2];
}

const double*
CCubicGrid::values() const
{
    return _values.data();
}

double*
CCubicGrid::values()
{
    return _values.data();
}
