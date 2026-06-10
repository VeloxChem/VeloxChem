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
#include <utility>

CSparseMatrix::CSparseMatrix()

    : _mtype(mat_t::general)

    , _block_rows{}

    , _block_cols{}

    , _blocks{}
{
}

CSparseMatrix::CSparseMatrix(const std::vector<CScreenedBasisFunctionPair> &screened_pairs, const mat_t mtype)

    : _mtype(mtype)

    , _block_rows{}

    , _block_cols{}

    , _blocks{}
{
    _block_rows.reserve(screened_pairs.size());

    _block_cols.reserve(screened_pairs.size());

    for (const auto &pair : screened_pairs)
    {
        const auto la = pair.bra_function().get_angular_momentum();

        const auto lb = pair.ket_function().get_angular_momentum();

        _block_rows.push_back((2 * la + 1) * (2 * lb + 1));

        _block_cols.push_back(static_cast<int>(pair.number_of_pairs()));
    }
}

CSparseMatrix::CSparseMatrix(const CSparseMatrix &other)

    : _mtype(other._mtype)

    , _block_rows(other._block_rows)

    , _block_cols(other._block_cols)

    , _blocks(other._blocks)
{
}

CSparseMatrix::CSparseMatrix(CSparseMatrix &&other) noexcept

    : _mtype(other._mtype)

    , _block_rows(std::move(other._block_rows))

    , _block_cols(std::move(other._block_cols))

    , _blocks(std::move(other._blocks))
{
}

auto
CSparseMatrix::operator=(const CSparseMatrix &other) -> CSparseMatrix &
{
    _mtype = other._mtype;

    _block_rows = other._block_rows;

    _block_cols = other._block_cols;

    _blocks = other._blocks;

    return *this;
}

auto
CSparseMatrix::operator=(CSparseMatrix &&other) noexcept -> CSparseMatrix &
{
    if (this != &other)
    {
        _mtype = other._mtype;

        _block_rows = std::move(other._block_rows);

        _block_cols = std::move(other._block_cols);

        _blocks = std::move(other._blocks);
    }

    return *this;
}

auto
CSparseMatrix::operator==(const CSparseMatrix &other) const -> bool
{
    return (_mtype == other._mtype) && (_block_rows == other._block_rows) && (_block_cols == other._block_cols) && (_blocks == other._blocks);
}

auto
CSparseMatrix::operator!=(const CSparseMatrix &other) const -> bool
{
    return !(*this == other);
}

auto
CSparseMatrix::type() const -> mat_t
{
    return _mtype;
}

auto
CSparseMatrix::number_of_keys() const -> size_t
{
    return _block_rows.size();
}

auto
CSparseMatrix::number_of_blocks() const -> size_t
{
    return _blocks.size();
}

auto
CSparseMatrix::has_block(const size_t key) const -> bool
{
    return _blocks.find(key) != _blocks.end();
}

auto
CSparseMatrix::block_rows(const size_t key) const -> int
{
    return _block_rows[key];
}

auto
CSparseMatrix::block_columns(const size_t key) const -> int
{
    return _block_cols[key];
}

auto
CSparseMatrix::block(const size_t key) -> CDenseMatrix &
{
    auto pos = _blocks.find(key);

    if (pos == _blocks.end())
    {
        CDenseMatrix mat(_block_rows[key], _block_cols[key]);

        mat.zero();

        pos = _blocks.emplace(key, std::move(mat)).first;
    }

    return pos->second;
}

auto
CSparseMatrix::block(const size_t key) const -> const CDenseMatrix &
{
    return _blocks.at(key);
}

auto
CSparseMatrix::keys() const -> std::vector<size_t>
{
    std::vector<size_t> block_keys;

    block_keys.reserve(_blocks.size());

    for (const auto &[key, mat] : _blocks) block_keys.push_back(key);

    std::ranges::sort(block_keys);

    return block_keys;
}

auto
CSparseMatrix::zero() -> void
{
    for (auto &[key, mat] : _blocks) mat.zero();
}
