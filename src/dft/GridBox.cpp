//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

CGridBox::CGridBox(const std::array<double, 6>& dimension, const CDenseMatrix& points)

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

auto
CGridBox::operator=(const CGridBox& source) -> CGridBox&
{
    if (this == &source) return *this;

    _dimension = source._dimension;

    _points = source._points;

    return *this;
}

auto
CGridBox::operator=(CGridBox&& source) noexcept -> CGridBox&
{
    if (this == &source) return *this;

    _dimension = std::move(source._dimension);

    _points = std::move(source._points);

    return *this;
}

auto
CGridBox::operator==(const CGridBox& other) const -> bool
{
    if (_dimension != other._dimension) return false;

    if (_points != other._points) return false;

    return true;
}

auto
CGridBox::operator!=(const CGridBox& other) const -> bool
{
    return !(*this == other);
}

auto
CGridBox::getBoxDimension() const -> std::array<double, 6>
{
    return _dimension;
}

auto
CGridBox::getGridPoints() const -> CDenseMatrix
{
    return _points;
}

auto
CGridBox::getNumberOfGridPoints() const -> int32_t
{
    return _points.getNumberOfColumns();
}

auto
CGridBox::getCoordinatesX() const -> const double*
{
    return _points.row(0);
}

auto
CGridBox::getCoordinatesY() const -> const double*
{
    return _points.row(1);
}

auto
CGridBox::getCoordinatesZ() const -> const double*
{
    return _points.row(2);
}

auto
CGridBox::getWeights() const -> const double*
{
    return _points.row(3);
}
