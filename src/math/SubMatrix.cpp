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

#include <algorithm>
#include <utility>

#include "CustomViews.hpp"
#include "MathFunc.hpp"

#include <iostream>

CSubMatrix::CSubMatrix()

    : _dimensions{0, 0, 0, 0}

    , _values{}
{
}

CSubMatrix::CSubMatrix(const std::array<size_t, 4> &dimensions)

    : _dimensions(dimensions)

    , _values(std::vector<double>(dimensions[2] * dimensions[3]))
{
}

CSubMatrix::CSubMatrix(const std::array<size_t, 4> &dimensions, const double factor)

    : _dimensions(dimensions)

    , _values(std::vector<double>(dimensions[2] * dimensions[3], factor))
{
}

CSubMatrix::CSubMatrix(const std::vector<double> &values, const std::array<size_t, 4> &dimensions)

    : _dimensions(dimensions)

    , _values(values)
{
}

CSubMatrix::CSubMatrix(const std::vector<double>&                    values,
                       const std::vector<std::pair<size_t, size_t>>& mask,
                       const std::array<size_t, 4>&                  dimensions)

    : _dimensions(dimensions)

    , _values(std::vector<double>(dimensions[2] * dimensions[3], 0.0))
{
    size_t idx = 0;
    
    for (const auto& index : mask)
    {
        if (index.first == index.second)
        {
            at(index) = values[idx];
        }
        else
        {
            at(index) = values[idx];
            
            at({index.second, index.first}) = values[idx];
        }
        
        idx++;
    }
}

CSubMatrix::CSubMatrix(const CSubMatrix &other)

    : _dimensions(other._dimensions)

    , _values(other._values)
{
}

CSubMatrix::CSubMatrix(CSubMatrix &&other) noexcept

    : _dimensions(std::move(other._dimensions))

    , _values(std::move(other._values))
{
}

auto
CSubMatrix::operator=(const CSubMatrix &other) -> CSubMatrix &
{
    _dimensions = other._dimensions;

    _values = other._values;

    return *this;
}

auto
CSubMatrix::operator=(CSubMatrix &&other) noexcept -> CSubMatrix &
{
    if (this != &other)
    {
        _dimensions = std::move(other._dimensions);

        _values = std::move(other._values);
    }

    return *this;
}

auto
CSubMatrix::operator==(const CSubMatrix &other) const -> bool
{
    if (_dimensions != other._dimensions)
    {
        return false;
    }
    else
    {
        return std::ranges::equal(_values, other._values, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); });
    }
}

auto
CSubMatrix::operator!=(const CSubMatrix &other) const -> bool
{
    return !(*this == other);
}

auto
CSubMatrix::operator+(const CSubMatrix &other) const -> CSubMatrix
{
    if (_dimensions == other._dimensions)
    {
        auto matrix = CSubMatrix(_dimensions);

        std::ranges::transform(_values, other._values, matrix._values.begin(), [](const double &a, const double &b) { return a + b; });

        return matrix;
    }
    else
    {
        return CSubMatrix();
    }
}

auto
CSubMatrix::set_offsets(const std::pair<size_t, size_t> &offsets) -> void
{
    _dimensions[0] = offsets.first;

    _dimensions[1] = offsets.second;
}

auto
CSubMatrix::set_values(const std::vector<double> &values) -> void
{
    _values = values;
}

auto
CSubMatrix::zero() -> void
{
    std::ranges::fill(_values, 0.0);
}

auto
CSubMatrix::scale(const double factor) -> void
{
    std::ranges::for_each(_values, [=](double &val) { val *= factor; });
}

auto
CSubMatrix::symmetrize() -> void
{
    if (is_square())
    {
        std::ranges::for_each(views::triangular(_dimensions[2]), [&](const auto &index) {
            const auto rindex = std::pair{index.second, index.first};
            const auto fval   = at(index) + at(rindex);
            at(index)         = fval;
            at(rindex)        = fval;
        });
    }
}

auto
CSubMatrix::get_dimensions() const -> std::array<size_t, 4>
{
    return _dimensions;
}

auto
CSubMatrix::get_values() const -> std::vector<double>
{
    return _values;
}
