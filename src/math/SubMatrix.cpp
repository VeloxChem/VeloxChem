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

#include "SubMatrix.hpp"

#include <cmath>
#include <cstring>

#include "AngularMomentum.hpp"

CSubMatrix::CSubMatrix()

    : _values(nullptr)

    , _dimensions({0, 0, 0, 0})
{
}

CSubMatrix::CSubMatrix(const T4Index& dimensions)

    : _values(nullptr)

    , _dimensions(dimensions)
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));
    }
}

CSubMatrix::CSubMatrix(const std::vector<double>& values, const T4Index& dimensions)

    : _values(nullptr)

    , _dimensions(dimensions)
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));

        std::memcpy(_values, values.data(), nelements * sizeof(double));
    }
}

CSubMatrix::CSubMatrix(const CSubMatrix& source)

    : _values(nullptr)

    , _dimensions(source._dimensions)
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));

        std::memcpy(_values, source._values, nelements * sizeof(double));
    }
}

CSubMatrix::CSubMatrix(CSubMatrix&& source) noexcept

    : _values(nullptr)

    , _dimensions(std::move(source._dimensions))
{
    _values = source._values;

    source._values = nullptr;
}

CSubMatrix::~CSubMatrix()
{
    if (_values != nullptr) std::free(_values);

    _values = nullptr;
}

CSubMatrix&
CSubMatrix::operator=(const CSubMatrix& source)
{
    if (this == &source) return *this;

    if (_values != nullptr) std::free(_values);

    _values = nullptr;

    _dimensions = source._dimensions;

    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));

        std::memcpy(_values, source._values, nelements * sizeof(double));
    }

    return *this;
}

CSubMatrix&
CSubMatrix::operator=(CSubMatrix&& source) noexcept
{
    if (this == &source) return *this;

    if (_values != nullptr) std::free(_values);

    _values = source._values;

    source._values = nullptr;

    _dimensions = std::move(source._dimensions);

    return *this;
}

auto
CSubMatrix::operator==(const CSubMatrix& other) const -> bool
{
    if (_dimensions != other._dimensions) return false;

    for (int64_t i = 0; i < _dimensions[2] * _dimensions[3]; i++)
    {
        if (std::fabs(_values[i] - other._values[i]) > 1.0e-13) return false;
    }

    return true;
}

auto
CSubMatrix::operator!=(const CSubMatrix& other) const -> bool
{
    return !(*this == other);
}

auto
CSubMatrix::setOffsets(const int64_t row_offset, const int64_t col_offset) -> void
{
    _dimensions[0] = row_offset;

    _dimensions[1] = col_offset;
}

auto
CSubMatrix::setValues(const std::vector<double>& values) -> void
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements == values.size())
    {
        std::memcpy(_values, values.data(), nelements * sizeof(double));
    }
}

auto
CSubMatrix::zero() -> void
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
#pragma omp simd
        for (int64_t i = 0; i < nelements; i++)
        {
            _values[i] = 0.0;
        }
    }
}

auto
CSubMatrix::symmetrize() -> void
{
    if (_dimensions[2] == _dimensions[3])
    {
        const auto nrows = _dimensions[2];

        for (int64_t i = 0; i < nrows; i++)
        {
            for (int64_t j = i; j < nrows; j++)
            {
                const auto ij_off = i * nrows + j;

                const auto ji_off = j * nrows + i;

                const auto fval = _values[ij_off] + _values[ji_off];

                _values[ij_off] = fval;

                _values[ji_off] = fval;
            }
        }
    }
}

auto
CSubMatrix::getDimensions() const -> T4Index
{
    return _dimensions;
}

auto
CSubMatrix::getValues() const -> std::vector<double>
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        std::vector<double> values(nelements, 0.0);

        std::memcpy(values.data(), _values, nelements * sizeof(double));

        return values;
    }
    else
    {
        return std::vector<double>();
    }
}

auto
CSubMatrix::getData() -> double*
{
    return _values;
}

auto
CSubMatrix::getData() const -> const double*
{
    return _values;
}
