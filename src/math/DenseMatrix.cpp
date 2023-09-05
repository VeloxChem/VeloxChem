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

#include "DenseMatrix.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <utility>

CDenseMatrix::CDenseMatrix()
{
}

CDenseMatrix::CDenseMatrix(const std::vector<double>& values, const int64_t nRows, const int64_t nColumns)

    : _nRows(static_cast<int32_t>(nRows))

    , _nColumns(static_cast<int32_t>(nColumns))

    , _values(values)
{
}

CDenseMatrix::CDenseMatrix(const int64_t nRows, const int64_t nColumns)

    : _nRows(static_cast<int32_t>(nRows))

    , _nColumns(static_cast<int32_t>(nColumns))

    , _values(std::vector<double>(static_cast<int32_t>(nRows * nColumns), 0.0))
{
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

CDenseMatrix&
CDenseMatrix::operator=(const CDenseMatrix& source)
{
    if (this == &source) return *this;

    _nRows = source._nRows;

    _nColumns = source._nColumns;

    _values = source._values;

    return *this;
}

CDenseMatrix&
CDenseMatrix::operator=(CDenseMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _nRows = std::move(source._nRows);

    _nColumns = std::move(source._nColumns);

    _values = std::move(source._values);

    return *this;
}

auto
CDenseMatrix::zero() -> void
{
    std::fill(_values.data(), _values.data() + _nRows * _nColumns, 0.0);
}

auto
CDenseMatrix::transpose() const -> CDenseMatrix
{
    CDenseMatrix tmat(_nColumns, _nRows);

    auto cvals = _values.data();

    auto tvals = tmat.values();

    for (int32_t i = 0; i < _nRows; i++)
    {
        for (int32_t j = 0; j < _nColumns; j++)
        {
            tvals[j * _nRows + i] = cvals[i * _nColumns + j];
        }
    }

    return tmat;
}

auto
CDenseMatrix::symmetrize() -> void
{
    auto fmat = _values.data();

    if (_nRows == _nColumns)
    {
        for (int32_t i = 0; i < _nRows; i++)
        {
            for (int32_t j = i; j < _nRows; j++)
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
    auto fmat = _values.data();

    if (_nRows == _nColumns)
    {
        for (int32_t i = 0; i < _nRows; i++)
        {
            for (int32_t j = i; j < _nRows; j++)
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
CDenseMatrix::getNumberOfRows() const -> int32_t
{
    return _nRows;
}

auto
CDenseMatrix::getNumberOfColumns() const -> int32_t
{
    return _nColumns;
}

auto
CDenseMatrix::getNumberOfElements() const -> int32_t
{
    return _nRows * _nColumns;
}

auto
CDenseMatrix::values() const -> const double*
{
    return _values.data();
}

auto
CDenseMatrix::values() -> double*
{
    return _values.data();
}

auto
CDenseMatrix::row(const int32_t iRow) const -> const double*
{
    if (iRow < getNumberOfRows())
    {
        return _values.data() + iRow * _nColumns;
    }
    else
    {
        return nullptr;
    }
}

auto
CDenseMatrix::row(const int32_t iRow) -> double*
{
    if (iRow < getNumberOfRows())
    {
        return _values.data() + iRow * _nColumns;
    }
    else
    {
        return nullptr;
    }
}
