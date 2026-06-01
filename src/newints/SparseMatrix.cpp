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
#include <map>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"

namespace newints {

namespace {

/// @brief Transposes a full block, optionally negating it (for antisymmetric
/// canonicalization of an i > j block to (j, i)).
auto
transpose_block(const Block &block, const bool negate) -> Block
{
    Block result{block.ncols, block.nrows, std::vector<double>(block.data.size()), Kind::full};

    const auto sign = negate ? -1.0 : 1.0;

    for (std::size_t r = 0; r < block.nrows; r++)
    {
        for (std::size_t c = 0; c < block.ncols; c++)
        {
            result.data[c * block.nrows + r] = sign * block.data[r * block.ncols + c];
        }
    }

    return result;
}

/// @brief Packs the lower triangle (r >= c) of a full square block into a
/// Kind::lower_triangular block of n(n+1)/2 values.
auto
pack_lower(const Block &block) -> Block
{
    const auto n = block.nrows;  // square diagonal block: nrows == ncols

    Block result{n, n, std::vector<double>(n * (n + 1) / 2), Kind::lower_triangular};

    for (std::size_t r = 0; r < n; r++)
    {
        for (std::size_t c = 0; c <= r; c++)
        {
            result.data[r * (r + 1) / 2 + c] = block.data[r * block.ncols + c];
        }
    }

    return result;
}

/// @brief Number of stored payload values for a block of the given dimensions
/// and layout (n(n+1)/2 packed for a lower-triangular diagonal block, else
/// nrows * ncols).
auto
payload_size(const std::size_t nrows, const std::size_t ncols, const Kind kind) -> std::size_t
{
    return (kind == Kind::lower_triangular) ? nrows * (nrows + 1) / 2 : nrows * ncols;
}

}  // namespace

auto
Block::operator==(const Block &other) const -> bool
{
    return (nrows == other.nrows) && (ncols == other.ncols) && (kind == other.kind) && (data == other.data);
}

SparseMatrix::SparseMatrix()

    : _symmetry(SymmetryType::general)

    , _meta{}

    , _sorted(true)

    , _data{}
{
}

SparseMatrix::SparseMatrix(const SymmetryType symmetry)

    : _symmetry(symmetry)

    , _meta{}

    , _sorted(true)

    , _data{}
{
}

auto
SparseMatrix::ensure_sorted() const -> void
{
    if (_sorted) return;

    std::ranges::stable_sort(_meta, [](const BlockMeta &a, const BlockMeta &b) { return (a.i != b.i) ? (a.i < b.i) : (a.j < b.j); });

    // collapse duplicate keys, keeping the most recently added (stable_sort kept
    // insertion order within an equal-key run, so the last occurrence wins)
    std::vector<BlockMeta> unique;

    unique.reserve(_meta.size());

    for (const auto &m : _meta)
    {
        if (!unique.empty() && unique.back().i == m.i && unique.back().j == m.j)
        {
            unique.back() = m;
        }
        else
        {
            unique.push_back(m);
        }
    }

    _meta.swap(unique);

    _sorted = true;
}

auto
SparseMatrix::operator==(const SparseMatrix &other) const -> bool
{
    if (_symmetry != other._symmetry) return false;

    ensure_sorted();

    other.ensure_sorted();

    if (_meta.size() != other._meta.size()) return false;

    for (std::size_t k = 0; k < _meta.size(); k++)
    {
        const auto &a = _meta[k];

        const auto &b = other._meta[k];

        if (a.i != b.i || a.j != b.j || a.nrows != b.nrows || a.ncols != b.ncols || a.kind != b.kind) return false;

        const auto n = payload_size(a.nrows, a.ncols, a.kind);

        if (!std::equal(_data.begin() + a.offset, _data.begin() + a.offset + n, other._data.begin() + b.offset)) return false;
    }

    return true;
}

auto
SparseMatrix::set_symmetry(const SymmetryType symmetry) -> void
{
    _symmetry = symmetry;
}

auto
SparseMatrix::add(const Key &key, const Block &block) -> void
{
    const auto [i, j] = key;

    // canonicalize to the stored block, then append its payload to the arena
    auto append = [this](const int bi, const int bj, const Block &b) {
        _meta.push_back(BlockMeta{bi, bj, static_cast<std::uint32_t>(b.nrows), static_cast<std::uint32_t>(b.ncols), _data.size(), b.kind});

        _data.insert(_data.end(), b.data.begin(), b.data.end());

        _sorted = false;
    };

    if (_symmetry == SymmetryType::general || i < j)
    {
        append(i, j, block);
    }
    else if (i > j)
    {
        append(j, i, transpose_block(block, _symmetry == SymmetryType::antisymmetric));
    }
    else  // i == j: diagonal block stored packed lower-triangular
    {
        append(i, j, pack_lower(block));
    }
}

auto
SparseMatrix::add(const int i, const int j, const Block &block) -> void
{
    add(Key{i, j}, block);
}

auto
SparseMatrix::reserve(const std::size_t count) -> void
{
    _meta.reserve(count);
}

auto
SparseMatrix::reserve_data(const std::size_t values) -> void
{
    _data.reserve(values);
}

auto
SparseMatrix::add_raw(const int i, const int j, const std::size_t nrows, const std::size_t ncols, const Kind kind, const double *src) -> void
{
    const auto n = payload_size(nrows, ncols, kind);

    _meta.push_back(BlockMeta{i, j, static_cast<std::uint32_t>(nrows), static_cast<std::uint32_t>(ncols), _data.size(), kind});

    _data.insert(_data.end(), src, src + n);

    _sorted = false;
}

auto
SparseMatrix::zero() -> void
{
    std::ranges::fill(_data, 0.0);
}

auto
SparseMatrix::symmetry() const -> SymmetryType
{
    return _symmetry;
}

auto
SparseMatrix::contains(const Key &key) const -> bool
{
    ensure_sorted();

    return std::ranges::binary_search(_meta, key, {}, [](const BlockMeta &m) { return Key{m.i, m.j}; });
}

auto
SparseMatrix::block(const Key &key) const -> std::optional<Block>
{
    ensure_sorted();

    const auto it = std::ranges::lower_bound(_meta, key, {}, [](const BlockMeta &m) { return Key{m.i, m.j}; });

    if (it == _meta.end() || it->i != key.first || it->j != key.second) return std::nullopt;

    const auto n = payload_size(it->nrows, it->ncols, it->kind);

    return Block{it->nrows, it->ncols, std::vector<double>(_data.begin() + it->offset, _data.begin() + it->offset + n), it->kind};
}

auto
SparseMatrix::number_of_blocks() const -> std::size_t
{
    ensure_sorted();

    return _meta.size();
}

auto
SparseMatrix::keys() const -> std::vector<Key>
{
    ensure_sorted();

    std::vector<Key> ckeys;

    ckeys.reserve(_meta.size());

    for (const auto &m : _meta) ckeys.push_back(Key{m.i, m.j});

    return ckeys;
}

auto
SparseMatrix::to_dense(const CMolecularBasis &basis) const -> CDenseMatrix
{
    const auto naos = static_cast<int>(basis.dimensions_of_basis());

    const auto maxl = basis.max_angular_momentum();

    // per-l number of contracted GTOs and starting offset in VeloxChem AO order
    std::vector<int> nbf(maxl + 1, 0);

    std::vector<int> base(maxl + 1, 0);

    for (int l = 0, acc = 0; l <= maxl; l++)
    {
        nbf[l]  = static_cast<int>(basis.number_of_basis_functions(l));
        base[l] = acc;
        acc += (2 * l + 1) * nbf[l];
    }

    // map each contracted-GTO (newints) index to the VeloxChem AO index of its components;
    // VeloxChem AO index = base[l] + m * nbf[l] + (#momentum-l shells on prior atoms + contraction)
    std::map<int, std::vector<int>> ao_map;

    const auto basis_indices = basis.basis_sets_indices();

    const auto atom_bases = basis.basis_sets();

    const auto natoms = static_cast<int>(basis_indices.size());

    std::vector<int> lcount_before(maxl + 1, 0);  // momentum-l shells on atoms processed so far

    int offset = 0;  // running newints (angular-expanded, atom-major) index

    for (int atom = 0; atom < natoms; atom++)
    {
        const auto &atom_basis = atom_bases[basis_indices[atom]];

        std::vector<int> c_within(maxl + 1, 0);  // contraction index within (atom, l)

        for (const auto &bf : atom_basis.basis_functions())
        {
            const auto l = bf.get_angular_momentum();

            const auto within = lcount_before[l] + c_within[l];

            std::vector<int> comps(2 * l + 1);

            for (int m = 0; m < 2 * l + 1; m++)
            {
                comps[m] = base[l] + m * nbf[l] + within;
            }

            ao_map[offset] = comps;

            c_within[l] += 1;

            offset += 2 * l + 1;
        }

        for (int l = 0; l <= maxl; l++) lcount_before[l] += c_within[l];
    }

    CDenseMatrix dmat(naos, naos);

    dmat.zero();

    const auto sign = (_symmetry == SymmetryType::antisymmetric) ? -1.0 : 1.0;

    ensure_sorted();

    for (const auto &m : _meta)
    {
        const auto &rcomps = ao_map.at(m.i);

        const auto &ccomps = ao_map.at(m.j);

        const auto *vals = _data.data() + m.offset;

        if (m.kind == Kind::lower_triangular)
        {
            // diagonal block (i == j): unpack the lower triangle, mirror with the symmetry sign
            const auto n = static_cast<int>(rcomps.size());

            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c <= r; c++)
                {
                    const auto v = vals[r * (r + 1) / 2 + c];

                    dmat.row(rcomps[r])[ccomps[c]] = v;

                    if (r != c) dmat.row(rcomps[c])[ccomps[r]] = sign * v;
                }
            }
        }
        else
        {
            const auto nr = static_cast<int>(rcomps.size());

            const auto nc = static_cast<int>(ccomps.size());

            for (int r = 0; r < nr; r++)
            {
                for (int c = 0; c < nc; c++)
                {
                    const auto v = vals[r * nc + c];

                    dmat.row(rcomps[r])[ccomps[c]] = v;

                    // for symmetric / antisymmetric matrices fill the implied partner (j, i)
                    if (_symmetry != SymmetryType::general) dmat.row(ccomps[c])[rcomps[r]] = sign * v;
                }
            }
        }
    }

    return dmat;
}

}  // namespace newints
