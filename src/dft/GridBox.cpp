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

#include "GridBox.hpp"

CGridBox::CGridBox()
{
}

CGridBox::CGridBox(const std::array<double, 6>& dimension, const CMemBlock2D<double>& points)

    : _dimension(dimension)

    , _points(points)
{
}

CGridBox::CGridBox(const CGridBox& source)

    : _dimension(source._dimension)

    , _points(source._points)
{
}

CGridBox::CGridBox(CGridBox&& source) noexcept

    : _dimension(std::move(source._dimension))

    , _points(std::move(source._points))
{
}

CGridBox::~CGridBox()
{
}

CGridBox&
CGridBox::operator=(const CGridBox& source)
{
    if (this == &source) return *this;

    _dimension = source._dimension;

    _points = source._points;

    return *this;
}

CGridBox&
CGridBox::operator=(CGridBox&& source) noexcept
{
    if (this == &source) return *this;

    _dimension = std::move(source._dimension);

    _points = std::move(source._points);

    return *this;
}

bool
CGridBox::operator==(const CGridBox& other) const
{
    if (_dimension != other._dimension) return false;

    if (_points != other._points) return false;

    return true;
}

bool
CGridBox::operator!=(const CGridBox& other) const
{
    return !(*this == other);
}

std::array<double, 6>
CGridBox::getBoxDimension() const
{
    return _dimension;
}

CMemBlock2D<double>
CGridBox::getGridPoints() const
{
    return _points;
}

int32_t
CGridBox::getNumberOfGridPoints() const
{
    return _points.size(0);
}

const double*
CGridBox::getCoordinatesX() const
{
    return _points.data(0);
}

const double*
CGridBox::getCoordinatesY() const
{
    return _points.data(1);
}

const double*
CGridBox::getCoordinatesZ() const
{
    return _points.data(2);
}

const double*
CGridBox::getWeights() const
{
    return _points.data(3);
}
