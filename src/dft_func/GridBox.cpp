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
CGridBox::getNumberOfGridPoints() const -> int
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
