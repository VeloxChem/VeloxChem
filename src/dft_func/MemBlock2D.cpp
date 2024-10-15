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

#include "MemBlock2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

CMemBlock2D::CMemBlock2D()

    : _nElements(0)
{
}

CMemBlock2D::CMemBlock2D(const int nElements, const int nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);

    _setDimensions();

    _data = std::vector<double>(_nElements);

    std::fill(_data.begin(), _data.end(), 0.0);
}

CMemBlock2D::CMemBlock2D(const CMemBlock2D& source)

    : _data(source._data)

    , _originalSizes(source._originalSizes)

    , _paddedSizes(source._paddedSizes)

    , _positions(source._positions)

    , _nElements(source._nElements)
{
}

CMemBlock2D::CMemBlock2D(CMemBlock2D&& source) noexcept

    : _data(std::move(source._data))

    , _originalSizes(std::move(source._originalSizes))

    , _paddedSizes(std::move(source._paddedSizes))

    , _positions(std::move(source._positions))

    , _nElements(std::move(source._nElements))
{
}

CMemBlock2D::~CMemBlock2D()
{
}

CMemBlock2D&
CMemBlock2D::operator=(const CMemBlock2D& source)
{
    if (this == &source) return *this;

    _nElements = source._nElements;

    _positions = source._positions;

    _paddedSizes = source._paddedSizes;

    _originalSizes = source._originalSizes;

    _data = source._data;

    return *this;
}

CMemBlock2D&
CMemBlock2D::operator=(CMemBlock2D&& source) noexcept
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
CMemBlock2D::operator==(const CMemBlock2D& other) const
{
    if (_nElements != other._nElements) return false;

    if (_positions != other._positions) return false;

    if (_paddedSizes != other._paddedSizes) return false;

    if (_originalSizes != other._originalSizes) return false;

    if (_data != other._data) return false;

    return true;
}

bool
CMemBlock2D::operator!=(const CMemBlock2D& other) const
{
    return !(*this == other);
}

void
CMemBlock2D::zero()
{
    std::fill(_data.begin(), _data.end(), 0.0);
}

bool
CMemBlock2D::isEmpty() const
{
    if (_nElements == 0) return true;

    return false;
}

double*
CMemBlock2D::data(const int iBlock)
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
CMemBlock2D::data(const int iBlock) const
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
CMemBlock2D::data(const int iBlock, const int iElement)
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
CMemBlock2D::data(const int iBlock, const int iElement) const
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
CMemBlock2D::size(const int iBlock) const
{
    if (iBlock < _originalSizes.size()) return _originalSizes.at(iBlock);

    return 0;
}

int
CMemBlock2D::pitched_size(const int iBlock) const
{
    if (iBlock < _paddedSizes.size()) return _paddedSizes.at(iBlock);

    return 0;
}

int
CMemBlock2D::blocks() const
{
    return _originalSizes.size();
}

void
CMemBlock2D::_setOriginalSizes(const int nElements, const int nBlocks)
{
    _originalSizes = std::vector<int>(nBlocks);

    for (int i = 0; i < nBlocks; i++)
        _originalSizes.at(i) = nElements;
}

void
CMemBlock2D::_setDimensions()
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
