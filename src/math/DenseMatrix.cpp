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

CDenseMatrix::CDenseMatrix(const std::vector<double>& values, const int32_t nRows, const int32_t nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)
{
    _values = std::vector<double>(nRows * nColumns);

    std::memcpy(_values.data(), values.data(), nRows * nColumns * sizeof(double));
}

CDenseMatrix::CDenseMatrix(const int32_t nRows, const int32_t nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)
{
    _values = std::vector<double>(nRows * nColumns);

    zero();
}

CDenseMatrix::CDenseMatrix(const CDenseMatrix& other)

    : _nRows(other.getNumberOfRows())

    , _nColumns(other.getNumberOfColumns())
{
    _values = std::vector<double>(other.getNumberOfElements());

    std::memcpy(_values.data(), other.values(), other.getNumberOfElements() * sizeof(double));
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
    std::fill(values(), values() + _nRows * _nColumns, 0.0);
}

auto
CDenseMatrix::transpose() const -> CDenseMatrix
{
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    CDenseMatrix tmat(ncols, nrows);

    auto cvals = values();

    auto tvals = tmat.values();

    for (int32_t i = 0; i < nrows; i++)
    {
        for (int32_t j = 0; j < ncols; j++)
        {
            tvals[j * nrows + i] = cvals[i * ncols + j];
        }
    }

    return tmat;
}

auto
CDenseMatrix::symmetrize() -> void
{
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    auto fmat = values();

    if (nrows == ncols)
    {
        for (int32_t i = 0; i < nrows; i++)
        {
            for (int32_t j = i; j < nrows; j++)
            {
                auto ijoff = i * ncols + j;

                auto jioff = j * ncols + i;

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
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    auto fmat = values();

    if (nrows == ncols)
    {
        for (int32_t i = 0; i < nrows; i++)
        {
            for (int32_t j = i; j < nrows; j++)
            {
                auto ijoff = i * ncols + j;

                auto jioff = j * ncols + i;

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
        return _values.data() + iRow * getNumberOfColumns();
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
        return _values.data() + iRow * getNumberOfColumns();
    }
    else
    {
        return nullptr;
    }
}
