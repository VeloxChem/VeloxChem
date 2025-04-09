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

CDenseMatrix::CDenseMatrix(const int nRows, const int nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)
{
    _values = std::vector<double>(_nRows * _nColumns);
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
    std::fill(_values.begin(), _values.end(), 0.0);
}

auto
CDenseMatrix::transpose() const -> CDenseMatrix
{
    CDenseMatrix tmat(_nColumns, _nRows);

    auto cvals = _values.data();

    auto tvals = tmat.values();

    for (int i = 0; i < _nRows; i++)
    {
        for (int j = 0; j < _nColumns; j++)
        {
            tvals[j * _nRows + i] = cvals[i * _nColumns + j];
        }
    }

    return tmat;
}

auto
CDenseMatrix::scale(const double factor) -> void
{
    for (int i = 0; i < _nRows; i++)
    {
        auto ptr_i = row(i);

        for (int j = 0; j < _nColumns; j++)
        {
            ptr_i[j] *= factor;
        }
    }
}

auto
CDenseMatrix::symmetrize() -> void
{
    auto fmat = _values.data();

    if (_nRows == _nColumns)
    {
        for (int i = 0; i < _nRows; i++)
        {
            for (int j = i; j < _nRows; j++)
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
        for (int i = 0; i < _nRows; i++)
        {
            for (int j = i; j < _nRows; j++)
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
CDenseMatrix::getNumberOfRows() const -> int
{
    return _nRows;
}

auto
CDenseMatrix::getNumberOfColumns() const -> int
{
    return _nColumns;
}

auto
CDenseMatrix::getNumberOfElements() const -> int
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
CDenseMatrix::row(const int iRow) const -> const double*
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
CDenseMatrix::row(const int iRow) -> double*
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
CDenseMatrix::slice(const int iPosition, const int nElements) const -> CDenseMatrix
{
    CDenseMatrix sliced_data(_nRows, nElements);

    for (int i = 0; i < _nRows; i++)
    {
        // set up pointers to data chunks

        auto idata = row(i) + iPosition;

        auto odata = sliced_data.row(i);

        // copy elements of data chunks

        for (int j = 0; j < nElements; j++)
        {
            odata[j] = idata[j];
        }
    }

    return sliced_data;
}
