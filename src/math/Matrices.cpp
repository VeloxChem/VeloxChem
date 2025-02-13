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

#include "Matrices.hpp"

#include <algorithm>

CMatrices::CMatrices()

    : _matrices(std::map<std::string, CMatrix*>())
{
}

CMatrices::CMatrices(const std::map<std::string, CMatrix>& matrices)

    : _matrices(std::map<std::string, CMatrix*>())
{
    std::ranges::for_each(matrices, [&](const auto& mvalue) { _matrices.insert({mvalue.first, new CMatrix(mvalue.second)}); });
}

CMatrices::CMatrices(const CMatrices& other)

    : _matrices(std::map<std::string, CMatrix*>())
{
    std::ranges::for_each(other._matrices, [&](const auto& mvalue) { _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)}); });
}

CMatrices::CMatrices(CMatrices&& other) noexcept

    : _matrices(std::move(other._matrices))
{
}

CMatrices::~CMatrices()
{
    _deallocate();
}

auto
CMatrices::operator=(const CMatrices& other) -> CMatrices&
{
    _deallocate();

    std::ranges::for_each(other._matrices, [&](const auto& mvalue) { _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)}); });

    return *this;
}

auto
CMatrices::operator=(CMatrices&& other) noexcept -> CMatrices&
{
    if (this != &other)
    {
        // release existing resources
        _deallocate();

        _matrices = std::move(other._matrices);
    }

    return *this;
}

auto
CMatrices::operator==(const CMatrices& other) const -> bool
{
    return std::ranges::equal(
        _matrices, other._matrices, [&](auto lhs, auto rhs) -> bool { return (lhs.first == rhs.first) && (*lhs.second == *rhs.second); });
}

auto
CMatrices::operator!=(const CMatrices& other) const -> bool
{
    return !(*this == other);
}

auto
CMatrices::add(const CMatrix& matrix, const std::string& key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const int key) -> void
{
    _matrices.insert({std::to_string(key), new CMatrix(matrix)});
}

auto
CMatrices::zero() -> void
{
    std::ranges::for_each(_matrices, [&](auto& mvalue) { mvalue.second->zero(); });
}

auto
CMatrices::scale(const double factor) -> void
{
    std::ranges::for_each(_matrices, [&](auto& mvalue) { mvalue.second->scale(factor); });
}

auto
CMatrices::symmetrize() -> void
{
    std::ranges::for_each(_matrices, [&](auto& mvalue) { mvalue.second->symmetrize(); });
}

auto
CMatrices::keys() const -> std::vector<std::string>
{
    std::vector<std::string> keys;

    keys.reserve(_matrices.size());

    std::ranges::for_each(_matrices, [&](auto& mvalue) { keys.push_back(mvalue.first); });

    return keys;
}

auto
CMatrices::matrix(const std::string& key) -> CMatrix*
{
    if (_matrices.contains(key))
    {
        return _matrices.at(key);
    }
    else
    {
        return nullptr;
    }
}

auto
CMatrices::matrix(const std::string& key) const -> const CMatrix*
{
    if (_matrices.contains(key))
    {
        return _matrices.at(key);
    }
    else
    {
        return nullptr;
    }
}

auto
CMatrices::matrix(const int key) -> CMatrix*
{
    if (const auto strkey = std::to_string(key); _matrices.contains(strkey))
    {
        return _matrices.at(strkey);
    }
    else
    {
        return nullptr;
    }
}

auto
CMatrices::matrix(const int key) const -> const CMatrix*
{
    if (const auto strkey = std::to_string(key); _matrices.contains(strkey))
    {
        return _matrices.at(strkey);
    }
    else
    {
        return nullptr;
    }
}

auto
CMatrices::_deallocate() -> void
{
    std::ranges::for_each(_matrices, [&](auto& mvalue) { delete mvalue.second; });

    _matrices.clear();
}
