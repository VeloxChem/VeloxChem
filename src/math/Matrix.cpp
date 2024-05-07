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

#include "Matrix.hpp"

#include <algorithm>
#include <cmath>

#include "AngularMomentum.hpp"

CMatrix::CMatrix()

    : _sub_matrices(std::map<T2Pair, CSubMatrix*>())

    , _mat_type(mat_t::gen)
{
}

CMatrix::CMatrix(const std::map<T2Pair, CSubMatrix>& sub_matrices, const mat_t mat_type)

    : _sub_matrices(std::map<T2Pair, CSubMatrix*>())

    , _mat_type(mat_type)
{
    for (const auto& mvalue : sub_matrices)
    {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(mvalue.second)});
    }
}

CMatrix::CMatrix(const CMatrix& other)

    : _sub_matrices(std::map<T2Pair, CSubMatrix*>())

    , _mat_type(other._mat_type)
{
    for (const auto& mvalue : other._sub_matrices)
    {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(*mvalue.second)});
    }
}

CMatrix::~CMatrix()
{
    for (auto& mvalue : _sub_matrices)
    {
        delete mvalue.second;
    }
}

auto
CMatrix::operator=(const CMatrix& source) -> CMatrix&
{
    if (this == &source) return *this;

    _mat_type = source._mat_type;

    for (auto& mvalue : _sub_matrices)
    {
        delete mvalue.second;
    }

    _sub_matrices.clear();

    for (const auto& mvalue : source._sub_matrices)
    {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(*mvalue.second)});
    }

    return *this;
}

auto
CMatrix::add(const CSubMatrix& sub_matrix, const T2Pair& angpair) -> void
{
    _sub_matrices.insert({angpair, new CSubMatrix(sub_matrix)});
}

auto
CMatrix::add(const T4Index& dimensions, const T2Pair& angpair) -> void
{
    _sub_matrices.insert({angpair, new CSubMatrix(dimensions)});
}

auto
CMatrix::setType(const mat_t mat_type) -> void
{
    _mat_type = mat_type;
}

auto
CMatrix::zero() -> void
{
    for (auto& mvalue : _sub_matrices)
    {
        mvalue.second->zero();
    }
}

auto
CMatrix::getAngularPairs() const -> std::vector<T2Pair>
{
    if (!_sub_matrices.empty())
    {
        std::vector<T2Pair> angpairs;

        for (const auto& mvalue : _sub_matrices)
        {
            angpairs.push_back(mvalue.first);
        }

        return angpairs;
    }
    else
    {
        return std::vector<T2Pair>();
    }
}

auto
CMatrix::getType() const -> mat_t
{
    return _mat_type;
}

auto
CMatrix::getSubMatrix(const T2Pair& angpair) -> CSubMatrix*
{
    for (const auto& mvalue : _sub_matrices)
    {
        if (mvalue.first == angpair)
        {
            return mvalue.second;
        }
    }

    if ((_mat_type == mat_t::symm) || (_mat_type == mat_t::antisymm))
    {
        const auto r_angpair = T2Pair({angpair.second, angpair.first});

        for (const auto& mvalue : _sub_matrices)
        {
            if (mvalue.first == r_angpair)
            {
                return mvalue.second;
            }
        }
    }

    return nullptr;
}

auto
CMatrix::getSubMatrix(const T2Pair& angpair) const -> const CSubMatrix*
{
    for (const auto& mvalue : _sub_matrices)
    {
        if (mvalue.first == angpair)
        {
            return mvalue.second;
        }
    }

    if ((_mat_type == mat_t::symm) || (_mat_type == mat_t::antisymm))
    {
        const auto r_angpair = T2Pair({angpair.second, angpair.first});

        for (const auto& mvalue : _sub_matrices)
        {
            if (mvalue.first == r_angpair)
            {
                return mvalue.second;
            }
        }
    }

    return nullptr;
}

auto
CMatrix::isAngularOrder(const T2Pair& angpair) const -> bool
{
    for (const auto& mvalue : _sub_matrices)
    {
        if (mvalue.first == angpair)
        {
            return true;
        }
    }

    return false;
}

auto
CMatrix::getNumberOfRows() const -> int64_t
{
    const auto row_keys = _getRowAngularKeys();

    const auto col_keys = _getColumnAngularKeys();

    if (row_keys.empty() || col_keys.empty())
    {
        return 0;
    }
    else
    {
        int64_t nrows = 0;

        const auto col_ang = *(col_keys.cbegin());

        for (const auto row_ang : row_keys)
        {
            if (const auto submat = getSubMatrix({row_ang, col_ang}); submat != nullptr)
            {
                if (isAngularOrder({row_ang, col_ang}))
                {
                    nrows += submat->getNumberOfRows();
                }
                else
                {
                    nrows += submat->getNumberOfColumns();
                }
            }
        }

        return nrows;
    }
}

auto
CMatrix::getNumberOfColumns() const -> int64_t
{
    const auto row_keys = _getRowAngularKeys();

    const auto col_keys = _getColumnAngularKeys();

    if (row_keys.empty() || col_keys.empty())
    {
        return 0;
    }
    else
    {
        int64_t ncols = 0;

        const auto row_ang = *(row_keys.cbegin());

        for (const auto col_ang : col_keys)
        {
            if (const auto submat = getSubMatrix({row_ang, col_ang}); submat != nullptr)
            {
                if (isAngularOrder({row_ang, col_ang}))
                {
                    ncols += submat->getNumberOfColumns();
                }
                else
                {
                    ncols += submat->getNumberOfRows();
                }
            }
        }

        return ncols;
    }
}

auto
CMatrix::getFullMatrix() const -> CSubMatrix
{
    const auto nrows = getNumberOfRows();

    const auto ncols = getNumberOfColumns();

    auto matrix = CSubMatrix({0, 0, nrows, ncols});

    for (const auto& mvalue : _sub_matrices)
    {
        const auto submat = mvalue.second;

        const auto [roff, coff, srows, scols] = submat->getDimensions();

        for (int64_t i = 0; i < srows; i++)
        {
            for (int64_t j = 0; j < scols; j++)
            {
                matrix.at(i + roff, j + coff, false) = submat->at(i, j, false);

                if (mvalue.first.first != mvalue.first.second)
                {
                    if (_mat_type == mat_t::symm)
                    {
                        matrix.at(j + coff, i + roff, false) = submat->at(i, j, false);
                    }

                    if (_mat_type == mat_t::antisymm)
                    {
                        matrix.at(j + coff, i + roff, false) = -submat->at(i, j, false);
                    }
                }
            }
        }
    }

    return matrix;
}

auto
CMatrix::_getRowAngularKeys() const -> std::set<int64_t>
{
    if (!_sub_matrices.empty())
    {
        std::set<int64_t> row_keys;

        for (const auto& mvalue : _sub_matrices)
        {
            row_keys.insert(mvalue.first.first);
        }

        return row_keys;
    }
    else
    {
        return std::set<int64_t>();
    }
}

auto
CMatrix::_getColumnAngularKeys() const -> std::set<int64_t>
{
    if (!_sub_matrices.empty())
    {
        std::set<int64_t> col_keys;

        for (const auto& mvalue : _sub_matrices)
        {
            col_keys.insert(mvalue.first.second);
        }

        return col_keys;
    }
    else
    {
        return std::set<int64_t>();
    }
}
