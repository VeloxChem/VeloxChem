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

const double*
CDense4DTensor::row(const int i, const int j, const int k) const
{
    if ((i < _iIndex) && (j < _jIndex) && (k < _kIndex))
    {
        return _values.data() + i * _jIndex * _kIndex * _lIndex + j * _kIndex * _lIndex + k * _lIndex;
    }
    else
    {
        return nullptr;
    }
}

double*
CDense4DTensor::row(const int i, const int j, const int k)
{
    if ((i < _iIndex) && (j < _jIndex) && (k < _kIndex))
    {
        return _values.data() + i * _jIndex * _kIndex * _lIndex + j * _kIndex * _lIndex + k * _lIndex;
    }
    else
    {
        return nullptr;
    }
}
