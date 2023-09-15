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
#include <sstream>
#include <utility>

CDenseMatrix::CDenseMatrix()
{
}

CDenseMatrix::CDenseMatrix(const std::vector<double>& values, const int32_t nRows, const int32_t nColumns)
{
    T4Index dim({0, 0, static_cast<int64_t>(nRows), static_cast<int64_t>(nColumns)});

    _values = CSubMatrix(values, dim);
}

CDenseMatrix::CDenseMatrix(const int32_t nRows, const int32_t nColumns)

{
    T4Index dim({0, 0, static_cast<int64_t>(nRows), static_cast<int64_t>(nColumns)});

    _values = CSubMatrix(dim);
}

CDenseMatrix::CDenseMatrix(const CDenseMatrix& other)

    : _values(other._values)
{
}

CDenseMatrix::~CDenseMatrix()
{
}

void
CDenseMatrix::zero()
{
    _values.zero();
}

CDenseMatrix
CDenseMatrix::transpose() const
{
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    CDenseMatrix tmat(ncols, nrows);

    auto cvals = _values.getData();

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

void
CDenseMatrix::symmetrize()
{
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    if (nrows == ncols)
    {
        auto fmat = _values.getData();

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

void
CDenseMatrix::symmetrizeAndScale(const double factor)
{
    auto nrows = getNumberOfRows();

    auto ncols = getNumberOfColumns();

    if (nrows == ncols)
    {
        auto fmat = _values.getData();

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

int32_t
CDenseMatrix::getNumberOfRows() const
{
    return static_cast<int32_t>(_values.getNumberOfRows());
}

int32_t
CDenseMatrix::getNumberOfColumns() const
{
    return static_cast<int32_t>(_values.getNumberOfColumns());
}

int32_t
CDenseMatrix::getNumberOfElements() const
{
    return getNumberOfRows() * getNumberOfColumns();
}

const double*
CDenseMatrix::values() const
{
    return _values.getData();
}

double*
CDenseMatrix::values()
{
    return _values.getData();
}

const double*
CDenseMatrix::row(const int32_t iRow) const
{
    if (iRow < getNumberOfRows())
    {
        return _values.getData() + iRow * getNumberOfColumns();
    }
    else
    {
        return nullptr;
    }
}

double*
CDenseMatrix::row(const int32_t iRow)
{
    if (iRow < getNumberOfRows())
    {
        return _values.getData() + iRow * getNumberOfColumns();
    }
    else
    {
        return nullptr;
    }
}
