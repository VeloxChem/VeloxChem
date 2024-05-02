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

#include "DenseMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <utility>

CDenseMatrix::CDenseMatrix()

    : _nRows(0)

    , _nColumns(0)
{
}

CDenseMatrix::CDenseMatrix(const CSubMatrix& submat)

    : _values(submat)
{
    const auto dim = _values.getDimensions();

    _nRows = dim[2];

    _nColumns = dim[3];
}

CDenseMatrix::CDenseMatrix(const int64_t nRows, const int64_t nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)
{
    T4Index dim({0, 0, _nRows, _nColumns});

    _values = CSubMatrix(dim);
}

CDenseMatrix::CDenseMatrix(const CDenseMatrix& source)

    : _nRows(source._nRows)

    , _nColumns(source._nColumns)

    , _values(source._values)
{
}

CDenseMatrix::CDenseMatrix(CDenseMatrix&& source) noexcept

    : _nRows(std::move(source._nRows))

    , _nColumns(std::move(source._nColumns))

    , _values(std::move(source._values))
{
}

CDenseMatrix::~CDenseMatrix()
{
}

auto
CDenseMatrix::operator=(const CDenseMatrix& source) -> CDenseMatrix&
{
    if (this == &source) return *this;

    _nRows = source._nRows;

    _nColumns = source._nColumns;

    _values = source._values;

    return *this;
}

auto
CDenseMatrix::operator=(CDenseMatrix&& source) noexcept -> CDenseMatrix&
{
    if (this == &source) return *this;

    _nRows = std::move(source._nRows);

    _nColumns = std::move(source._nColumns);

    _values = std::move(source._values);

    return *this;
}

auto
CDenseMatrix::operator==(const CDenseMatrix& other) const -> bool
{
    if (_nRows != other._nRows) return false;

    if (_nColumns != other._nColumns) return false;

    if (_values != other._values) return false;

    return true;
}

auto
CDenseMatrix::operator!=(const CDenseMatrix& other) const -> bool
{
    return !(*this == other);
}

auto
CDenseMatrix::zero() -> void
{
    _values.zero();
}

auto
CDenseMatrix::transpose() const -> CDenseMatrix
{
    CDenseMatrix tmat(_nColumns, _nRows);

    auto cvals = _values.getData();

    auto tvals = tmat.values();

    for (int64_t i = 0; i < _nRows; i++)
    {
        for (int64_t j = 0; j < _nColumns; j++)
        {
            tvals[j * _nRows + i] = cvals[i * _nColumns + j];
        }
    }

    return tmat;
}

auto
CDenseMatrix::symmetrize() -> void
{
    auto fmat = _values.getData();

    if (_nRows == _nColumns)
    {
        for (int64_t i = 0; i < _nRows; i++)
        {
            for (int64_t j = i; j < _nRows; j++)
            {
                auto ijoff = i * _nColumns + j;

                auto jioff = j * _nColumns + i;

                auto fval = fmat[ijoff] + fmat[jioff];

                fmat[ijoff] = fval;

                fmat[jioff] = fval;
            }
        }
    }
}

auto
CDenseMatrix::symmetrizeAndScale(const double factor) -> void
{
    auto fmat = _values.getData();

    if (_nRows == _nColumns)
    {
        for (int64_t i = 0; i < _nRows; i++)
        {
            for (int64_t j = i; j < _nRows; j++)
            {
                auto ijoff = i * _nColumns + j;

                auto jioff = j * _nColumns + i;

                auto fval = factor * (fmat[ijoff] + fmat[jioff]);

                fmat[ijoff] = fval;

                fmat[jioff] = fval;
            }
        }
    }
}

auto
CDenseMatrix::getNumberOfRows() const -> int64_t
{
    return _nRows;
}

auto
CDenseMatrix::getNumberOfColumns() const -> int64_t
{
    return _nColumns;
}

auto
CDenseMatrix::getNumberOfElements() const -> int64_t
{
    return _nRows * _nColumns;
}

auto
CDenseMatrix::values() const -> const double*
{
    return _values.getData();
}

auto
CDenseMatrix::values() -> double*
{
    return _values.getData();
}

auto
CDenseMatrix::row(const int64_t iRow) const -> const double*
{
    if (iRow < getNumberOfRows())
    {
        return _values.getData() + iRow * _nColumns;
    }
    else
    {
        return nullptr;
    }
}

auto
CDenseMatrix::row(const int64_t iRow) -> double*
{
    if (iRow < getNumberOfRows())
    {
        return _values.getData() + iRow * _nColumns;
    }
    else
    {
        return nullptr;
    }
}

auto
CDenseMatrix::slice(const int64_t iPosition, const int64_t nElements) const -> CDenseMatrix
{
    CDenseMatrix sliced_data(_nRows, nElements);

    for (int64_t i = 0; i < _nRows; i++)
    {
        // set up pointers to data chunks

        auto idata = row(i) + iPosition;

        auto odata = sliced_data.row(i);

        // copy elements of data chunks

        for (int64_t j = 0; j < nElements; j++)
        {
            odata[j] = idata[j];
        }
    }

    return sliced_data;
}
