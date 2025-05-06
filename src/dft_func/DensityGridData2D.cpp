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

#include "DensityGridData2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

CDensityGridData2D::CDensityGridData2D()

    : _nElements(0)
{
}

CDensityGridData2D::CDensityGridData2D(const int nElements, const int nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);

    _setDimensions();

    _data = std::vector<double>(_nElements);

    std::fill(_data.begin(), _data.end(), 0.0);
}

CDensityGridData2D::CDensityGridData2D(const CDensityGridData2D& source)

    : _data(source._data)

    , _originalSizes(source._originalSizes)

    , _paddedSizes(source._paddedSizes)

    , _positions(source._positions)

    , _nElements(source._nElements)
{
}

CDensityGridData2D::CDensityGridData2D(CDensityGridData2D&& source) noexcept

    : _data(std::move(source._data))

    , _originalSizes(std::move(source._originalSizes))

    , _paddedSizes(std::move(source._paddedSizes))

    , _positions(std::move(source._positions))

    , _nElements(std::move(source._nElements))
{
}

CDensityGridData2D::~CDensityGridData2D()
{
}

CDensityGridData2D&
CDensityGridData2D::operator=(const CDensityGridData2D& source)
{
    if (this == &source) return *this;

    _nElements = source._nElements;

    _positions = source._positions;

    _paddedSizes = source._paddedSizes;

    _originalSizes = source._originalSizes;

    _data = source._data;

    return *this;
}

CDensityGridData2D&
CDensityGridData2D::operator=(CDensityGridData2D&& source) noexcept
{
    if (this == &source) return *this;

    _nElements = std::move(source._nElements);

    _positions = std::move(source._positions);

    _paddedSizes = std::move(source._paddedSizes);

    _originalSizes = std::move(source._originalSizes);

    _data = std::move(source._data);

    return *this;
}

bool
CDensityGridData2D::operator==(const CDensityGridData2D& other) const
{
    if (_nElements != other._nElements) return false;

    if (_positions != other._positions) return false;

    if (_paddedSizes != other._paddedSizes) return false;

    if (_originalSizes != other._originalSizes) return false;

    if (_data != other._data) return false;

    return true;
}

bool
CDensityGridData2D::operator!=(const CDensityGridData2D& other) const
{
    return !(*this == other);
}

void
CDensityGridData2D::zero()
{
    std::fill(_data.begin(), _data.end(), 0.0);
}

double*
CDensityGridData2D::data(const int iBlock)
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;

        if (iBlock >= blocks()) return nullptr;

        return _data.data() + _positions.at(iBlock);
    }

    return nullptr;
}

const double*
CDensityGridData2D::data(const int iBlock) const
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;

        if (iBlock >= blocks()) return nullptr;

        return _data.data() + _positions.at(iBlock);
    }

    return nullptr;
}

double*
CDensityGridData2D::data(const int iBlock, const int iElement)
{
    if (_originalSizes.size() > 0)
    {
        auto pdata = _data.data() + _positions.at(iBlock);

        if (iElement < _originalSizes.at(iBlock))
        {
            return &(pdata[iElement]);
        }

        return nullptr;
    }

    return nullptr;
}

const double*
CDensityGridData2D::data(const int iBlock, const int iElement) const
{
    if (_originalSizes.size() > 0)
    {
        auto pdata = _data.data() + _positions.at(iBlock);

        if (iElement < _originalSizes.at(iBlock))
        {
            return &(pdata[iElement]);
        }

        return nullptr;
    }

    return nullptr;
}

int
CDensityGridData2D::size(const int iBlock) const
{
    if (iBlock < _originalSizes.size()) return _originalSizes.at(iBlock);

    return 0;
}

int
CDensityGridData2D::pitched_size(const int iBlock) const
{
    if (iBlock < _paddedSizes.size()) return _paddedSizes.at(iBlock);

    return 0;
}

int
CDensityGridData2D::blocks() const
{
    return _originalSizes.size();
}

void
CDensityGridData2D::_setOriginalSizes(const int nElements, const int nBlocks)
{
    _originalSizes = std::vector<int>(nBlocks);

    for (int i = 0; i < nBlocks; i++)
        _originalSizes.at(i) = nElements;
}

void
CDensityGridData2D::_setDimensions()
{
    auto numblocks = _originalSizes.size();

    _paddedSizes = std::vector<int>(numblocks);

    _positions = std::vector<int>(numblocks);

    // loop over data chunks

    int primdim = 64 / sizeof(double);

    _nElements = 0;

    for (int i = 0; i < numblocks; i++)
    {
        // compute padded size of data chunk

        auto pblocks = _originalSizes.at(i) / primdim;

        if ((_originalSizes.at(i) % primdim) != 0) pblocks++;

        _paddedSizes.at(i) = pblocks * primdim;

        // determine start position of data chunk

        _positions.at(i) = _nElements;

        _nElements += _paddedSizes.at(i);
    }
}
