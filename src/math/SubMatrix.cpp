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

#include "SubMatrix.hpp"

#include <algorithm>
#include <utility>

#include "CustomViews.hpp"
#include "MathFunc.hpp"

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

CSubMatrix::CSubMatrix(const CSubMatrix &other)

    : _dimensions(other._dimensions)

    , _values(other._values)
{
}

CSubMatrix::CSubMatrix(CSubMatrix &&other) noexcept

    : _dimensions{0, 0, 0, 0}

    , _values{}
{
    std::swap(_dimensions, other._dimensions);

    std::swap(_values, other._values);
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
    std::swap(_dimensions, other._dimensions);

    std::swap(_values, other._values);

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
