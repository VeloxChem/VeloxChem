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
