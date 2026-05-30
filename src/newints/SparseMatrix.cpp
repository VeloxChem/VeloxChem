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

#include "SparseMatrix.hpp"

#include <algorithm>

namespace newints {

auto
Block::operator==(const Block &other) const -> bool
{
    return (nrows == other.nrows) && (ncols == other.ncols) && (data == other.data);
}

SparseMatrix::SparseMatrix()

    : _symmetry(SymmetryType::general)

    , _blocks{}
{
}

SparseMatrix::SparseMatrix(const SymmetryType symmetry)

    : _symmetry(symmetry)

    , _blocks{}
{
}

auto
SparseMatrix::operator==(const SparseMatrix &other) const -> bool
{
    return (_symmetry == other._symmetry) && (_blocks == other._blocks);
}

auto
SparseMatrix::set_symmetry(const SymmetryType symmetry) -> void
{
    _symmetry = symmetry;
}

auto
SparseMatrix::add(const Key &key, const Block &block) -> void
{
    _blocks[key] = block;
}

auto
SparseMatrix::add(const int i, const int j, const Block &block) -> void
{
    _blocks[{i, j}] = block;
}

auto
SparseMatrix::zero() -> void
{
    std::ranges::for_each(_blocks, [](auto &entry) { std::ranges::fill(entry.second.data, 0.0); });
}

auto
SparseMatrix::symmetry() const -> SymmetryType
{
    return _symmetry;
}

auto
SparseMatrix::contains(const Key &key) const -> bool
{
    return _blocks.contains(key);
}

auto
SparseMatrix::block(const Key &key) -> Block *
{
    auto it = _blocks.find(key);

    return (it != _blocks.end()) ? &it->second : nullptr;
}

auto
SparseMatrix::block(const Key &key) const -> const Block *
{
    auto it = _blocks.find(key);

    return (it != _blocks.end()) ? &it->second : nullptr;
}

auto
SparseMatrix::number_of_blocks() const -> std::size_t
{
    return _blocks.size();
}

auto
SparseMatrix::keys() const -> std::vector<Key>
{
    std::vector<Key> ckeys;

    ckeys.reserve(_blocks.size());

    std::ranges::for_each(_blocks, [&](const auto &entry) { ckeys.push_back(entry.first); });

    return ckeys;
}

}  // namespace newints
