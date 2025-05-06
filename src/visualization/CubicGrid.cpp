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

#include "CubicGrid.hpp"

#include <cstring>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"

CCubicGrid::CCubicGrid()

    : _origin{0.0, 0.0, 0.0}

    , _stepSize{0.0, 0.0, 0.0}

    , _numPoints{0, 0, 0}

    , _values{std::vector<double>()}
{
}

CCubicGrid::CCubicGrid(const std::array<double, 3>& origin, const std::array<double, 3>& stepSize, const std::array<int, 3>& numPoints)

    : _origin(origin)

    , _stepSize(stepSize)

    , _numPoints(numPoints)
{
    std::string errmsg("CubicGrid: Incorrect dimension");

    errors::assertMsgCritical(_origin.size() == 3, errmsg);

    errors::assertMsgCritical(_stepSize.size() == 3, errmsg);

    errors::assertMsgCritical(_numPoints.size() == 3, errmsg);

    _values = std::vector<double>(numPoints[0] * numPoints[1] * numPoints[2]);
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

std::array<double, 3>
CCubicGrid::getOrigin() const
{
    return _origin;
}

std::array<double, 3>
CCubicGrid::getStepSize() const
{
    return _stepSize;
}

std::array<int, 3>
CCubicGrid::getNumPoints() const
{
    return _numPoints;
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

void
CCubicGrid::setValues(const std::vector<double> vals)
{
    std::string errmsg("CubicGrid.set_values: Inconsistent number of grid points");

    auto npoints = _numPoints[0] * _numPoints[1] * _numPoints[2];

    errors::assertMsgCritical(npoints == static_cast<int>(vals.size()), errmsg);

    std::memcpy(_values.data(), vals.data(), npoints * sizeof(double));
}
