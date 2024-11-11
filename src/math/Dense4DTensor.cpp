//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "Dense4DTensor.hpp"

#include <cmath>
#include <utility>

#include "StringFormat.hpp"

CDense4DTensor::CDense4DTensor()

    : _iIndex(0)

    , _jIndex(0)

    , _kIndex(0)

    , _lIndex(0)
{
}

CDense4DTensor::CDense4DTensor(const std::vector<double>& values, const int iIndex, const int jIndex, const int kIndex, const int lIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)

    , _lIndex(lIndex)

    , _values(values)
{
}

CDense4DTensor::CDense4DTensor(const int iIndex, const int jIndex, const int kIndex, const int lIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)

    , _lIndex(lIndex)

    , _values(std::vector<double>(iIndex * jIndex * kIndex * lIndex, 0.0))
{
}

CDense4DTensor::CDense4DTensor(const int nRows)

    : _iIndex(nRows)

    , _jIndex(nRows)

    , _kIndex(nRows)

    , _lIndex(nRows)

    , _values(std::vector<double>(nRows * nRows * nRows * nRows, 0.0))
{
}

CDense4DTensor::CDense4DTensor(const CDense4DTensor& source)

    : _iIndex(source._iIndex)

    , _jIndex(source._jIndex)

    , _kIndex(source._lIndex)

    , _lIndex(source._lIndex)

    , _values(source._values)
{
}

CDense4DTensor::CDense4DTensor(CDense4DTensor&& source) noexcept

    : _iIndex(std::move(source._iIndex))

    , _jIndex(std::move(source._jIndex))

    , _kIndex(std::move(source._kIndex))

    , _lIndex(std::move(source._lIndex))

    , _values(std::move(source._values))
{
}

CDense4DTensor::~CDense4DTensor()
{
}

CDense4DTensor&
CDense4DTensor::operator=(const CDense4DTensor& source)
{
    if (this == &source) return *this;

    _iIndex = source._iIndex;

    _jIndex = source._jIndex;

    _lIndex = source._kIndex;

    _lIndex = source._lIndex;

    _values = source._values;

    return *this;
}

CDense4DTensor&
CDense4DTensor::operator=(CDense4DTensor&& source) noexcept
{
    if (this == &source) return *this;

    _iIndex = std::move(source._iIndex);

    _jIndex = std::move(source._jIndex);

    _kIndex = std::move(source._kIndex);

    _lIndex = std::move(source._lIndex);

    _values = std::move(source._values);

    return *this;
}

bool
CDense4DTensor::operator==(const CDense4DTensor& other) const
{
    if (_iIndex != other._iIndex) return false;

    if (_jIndex != other._jIndex) return false;

    if (_kIndex != other._kIndex) return false;

    if (_lIndex != other._lIndex) return false;

    if (_values != other._values) return false;

    return true;
}

bool
CDense4DTensor::operator!=(const CDense4DTensor& other) const
{
    return !(*this == other);
}

void
CDense4DTensor::zero()
{
    std::fill(_values.begin(), _values.end(), 0.0);
}

int
CDense4DTensor::getiIndex() const
{
    return _iIndex;
}

int
CDense4DTensor::getjIndex() const
{
    return _jIndex;
}

int
CDense4DTensor::getkIndex() const
{
    return _kIndex;
}

int
CDense4DTensor::getlIndex() const
{
    return _lIndex;
}

int
CDense4DTensor::getNumberOfElements() const
{
    return _iIndex * _jIndex * _kIndex * _lIndex;
}

const double*
CDense4DTensor::values() const
{
    return _values.data();
}

double*
CDense4DTensor::values()
{
    return _values.data();
}
