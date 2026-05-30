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

}  // namespace

auto
Block::operator==(const Block &other) const -> bool
{
    return (nrows == other.nrows) && (ncols == other.ncols) && (kind == other.kind) && (data == other.data);
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
    const auto [i, j] = key;

    if (_symmetry == SymmetryType::general)
    {
        _blocks[key] = block;
    }
    else if (i < j)
    {
        _blocks[key] = block;
    }
    else if (i > j)
    {
        _blocks[{j, i}] = transpose_block(block, _symmetry == SymmetryType::antisymmetric);
    }
    else  // i == j: diagonal block stored packed lower-triangular
    {
        _blocks[key] = pack_lower(block);
    }
}

auto
SparseMatrix::add(const int i, const int j, const Block &block) -> void
{
    add(Key{i, j}, block);
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

    for (const auto &[key, block] : _blocks)
    {
        const auto &rcomps = ao_map.at(key.first);

        const auto &ccomps = ao_map.at(key.second);

        if (block.kind == Kind::lower_triangular)
        {
            // diagonal block (i == j): unpack the lower triangle, mirror with the symmetry sign
            const auto n = static_cast<int>(rcomps.size());

            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c <= r; c++)
                {
                    const auto v = block.data[r * (r + 1) / 2 + c];

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
                    const auto v = block.data[r * nc + c];

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
