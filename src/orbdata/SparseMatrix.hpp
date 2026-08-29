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


#ifndef SparseMatrix_hpp
#define SparseMatrix_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairGroup.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/// @brief Defines supported states of the values blocks of a sparse matrix:
/// valstat::empty     - the values blocks are not allocated
/// valstat::allocated - the values blocks are allocated, with undefined content
enum class valstat
{
    empty,
    allocated
};

/// @brief Class CSparseMatrix stores the sparsity pattern of a matrix in atomic
/// orbital basis, as the sparsity patterns of its off-diagonal and diagonal atom
/// pair blocks. The off-diagonal blocks are described by CAtomBasisPairSparsity
/// and the diagonal blocks by CAtomBasisDiagonalSparsity, whose storage
/// requirements differ.
/// @note The blocks refer to their atom bases by index among the unique atom
/// bases of the molecular basis they originate from. The molecular basis is not
/// stored here, so which molecular basis an index refers to follows from the
/// matrix type: a symmetric or antisymmetric matrix has a single molecular basis
/// on both sides, while a general matrix has a molecular basis on each side,
/// with the bra index referring to the bra molecular basis and the ket index to
/// the ket molecular basis.
/// @note The blocks without atom pairs or atoms are dropped, so a block index is
/// not an index of the atom basis pair group it originates from, and the numbers
/// of off-diagonal and diagonal blocks need not agree.
/// @note The values blocks are not allocated by the constructors, which set up
/// the sparsity patterns only. The values blocks are addressed as one pointer per
/// block, with the off-diagonal blocks preceding the diagonal blocks, and the
/// state of the values blocks is reported by get_values_state().
class CSparseMatrix
{
   public:
    /// @brief The default constructor.
    CSparseMatrix()

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_t::general)

        , _values_state(valstat::empty)
    {
    }

    /// @brief The constructor with sparsity patterns of the blocks and matrix type.
    /// @param pair_blocks The sparsity patterns of the off-diagonal blocks.
    /// @param diagonal_blocks The sparsity patterns of the diagonal blocks.
    /// @param mat_type The type of matrix.
    CSparseMatrix(const std::vector<CAtomBasisPairSparsity>     &pair_blocks,
                  const std::vector<CAtomBasisDiagonalSparsity> &diagonal_blocks,
                  const mat_t                                    mat_type)

        : _pair_blocks(pair_blocks)

        , _diagonal_blocks(diagonal_blocks)

        , _values(pair_blocks.size() + diagonal_blocks.size(), nullptr)

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
    }

    /// @brief The constructor with molecule, molecular basis, integral screener
    /// and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param screener The integral bound, evaluated as screener(bra_function,
    /// ket_function, distance).
    /// @param threshold The screening threshold.
    /// @param mat_type The type of matrix.
    /// @param storage The storage layout of the diagonal blocks.
    template <typename B>
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage)

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        // NOTE: a symmetric or antisymmetric matrix needs the upper triangle of
        // the atom basis pair groups only, while a general matrix needs their
        // full direct product, which the two molecular bases factory delivers
        // even when handed the same molecular basis twice.

        auto groups = (mat_type == mat_t::general) ? basis.basis_pair_groups(basis) : basis.basis_pair_groups();

        _add_blocks(molecule, groups, screener, threshold, storage);

        _values.assign(_pair_blocks.size() + _diagonal_blocks.size(), nullptr);
    }

    /// @brief The constructor with molecule, molecular bases on bra and ket
    /// sides, integral screener and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param screener The integral bound, evaluated as screener(bra_function,
    /// ket_function, distance).
    /// @param threshold The screening threshold.
    /// @param mat_type The type of matrix, which must be general.
    /// @param storage The storage layout of the diagonal blocks.
    template <typename B>
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &bra_basis,
                  const CMolecularBasis &ket_basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage)

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        errors::assertMsgCritical(mat_type == mat_t::general,
                                  std::string("SparseMatrix: Matrix with two molecular bases must be of general type"));

        auto groups = bra_basis.basis_pair_groups(ket_basis);

        _add_blocks(molecule, groups, screener, threshold, storage);

        _values.assign(_pair_blocks.size() + _diagonal_blocks.size(), nullptr);
    }

    /// @brief The copy constructor.
    /// @param other The sparse matrix to be copied.
    CSparseMatrix(const CSparseMatrix &other)

        : _pair_blocks(other._pair_blocks)

        , _diagonal_blocks(other._diagonal_blocks)

        , _values(other._values.size(), nullptr)

        , _type(other._type)

        , _values_state(other._values_state)
    {
        if (_values_state == valstat::allocated) _copy_values(other);
    }

    /// @brief The move constructor.
    /// @param other The sparse matrix to be moved.
    CSparseMatrix(CSparseMatrix &&other) noexcept

        : _pair_blocks(std::move(other._pair_blocks))

        , _diagonal_blocks(std::move(other._diagonal_blocks))

        , _values(std::move(other._values))

        , _type(other._type)

        , _values_state(other._values_state)
    {
        other._values.clear();

        other._values_state = valstat::empty;
    }

    /// @brief The destructor.
    ~CSparseMatrix()
    {
        _deallocate();
    }

    /// @brief The copy assignment operator.
    /// @param other The sparse matrix to be copy assigned.
    /// @return The assigned sparse matrix.
    auto
    operator=(const CSparseMatrix &other) -> CSparseMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _pair_blocks = other._pair_blocks;

            _diagonal_blocks = other._diagonal_blocks;

            _values.assign(other._values.size(), nullptr);

            _type = other._type;

            _values_state = other._values_state;

            if (_values_state == valstat::allocated) _copy_values(other);
        }

        return *this;
    }

    /// @brief The move assignment operator.
    /// @param other The sparse matrix to be move assigned.
    /// @return The assigned sparse matrix.
    auto
    operator=(CSparseMatrix &&other) noexcept -> CSparseMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _pair_blocks = std::move(other._pair_blocks);

            _diagonal_blocks = std::move(other._diagonal_blocks);

            _values = std::move(other._values);

            _type = other._type;

            _values_state = other._values_state;

            other._values.clear();

            other._values_state = valstat::empty;
        }

        return *this;
    }

    /// @brief Allocates the values blocks of matrix, leaving their content
    /// undefined. Each values block is allocated separately, so no single
    /// allocation spans the whole matrix.
    auto
    allocate() -> void
    {
        if (_values_state == valstat::allocated) return;

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()),
                              [&](const auto i) { _values[i] = new double[number_of_elements(i)]; });

        _values_state = valstat::allocated;
    }

    /// @brief Sets the values of all values blocks of matrix to zero.
    auto
    zero() -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseMatrix.zero: Values blocks of matrix are not allocated"));

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()),
                              [&](const auto i) { std::fill(_values[i], _values[i] + number_of_elements(i), 0.0); });
    }

    /// @brief Deallocates the values blocks of matrix, leaving its sparsity
    /// patterns intact.
    auto
    deallocate() -> void
    {
        _deallocate();
    }

    /// @brief Scales the values of all values blocks of matrix.
    /// @param factor The scaling factor.
    auto
    scale(const double factor) -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseMatrix.scale: Values blocks of matrix are not allocated"));

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()), [&](const auto i) {
            std::ranges::for_each(std::views::iota(size_t{0}, number_of_elements(i)), [&](const auto k) { _values[i][k] *= factor; });
        });
    }

    /// @brief Gets values of the off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @return The pointer to the values of the block.
    auto
    pair_values(const size_t index) -> double *
    {
        _check_values(index, _pair_blocks.size(), std::string("pair_values"));

        return _values[index];
    }

    /// @brief Gets values of the off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @return The constant pointer to the values of the block.
    auto
    pair_values(const size_t index) const -> const double *
    {
        _check_values(index, _pair_blocks.size(), std::string("pair_values"));

        return _values[index];
    }

    /// @brief Gets values of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The pointer to the values of the block.
    auto
    diagonal_values(const size_t index) -> double *
    {
        _check_values(index, _diagonal_blocks.size(), std::string("diagonal_values"));

        return _values[_pair_blocks.size() + index];
    }

    /// @brief Gets values of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The constant pointer to the values of the block.
    auto
    diagonal_values(const size_t index) const -> const double *
    {
        _check_values(index, _diagonal_blocks.size(), std::string("diagonal_values"));

        return _values[_pair_blocks.size() + index];
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The pointer to the values of the combination of basis functions.
    auto
    pair_values(const size_t index,
                const int    bra_angular_momentum,
                const size_t bra_index,
                const int    ket_angular_momentum,
                const size_t ket_index) -> double *
    {
        return pair_values(index) +
               _pair_blocks[index].element_offset(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The constant pointer to the values of the combination of basis functions.
    auto
    pair_values(const size_t index,
                const int    bra_angular_momentum,
                const size_t bra_index,
                const int    ket_angular_momentum,
                const size_t ket_index) const -> const double *
    {
        return pair_values(index) +
               _pair_blocks[index].element_offset(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The pointer to the values of the combination of basis functions.
    auto
    diagonal_values(const size_t index,
                    const int    bra_angular_momentum,
                    const size_t bra_index,
                    const int    ket_angular_momentum,
                    const size_t ket_index) -> double *
    {
        return diagonal_values(index) +
               _diagonal_blocks[index].element_offset(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The constant pointer to the values of the combination of basis functions.
    auto
    diagonal_values(const size_t index,
                    const int    bra_angular_momentum,
                    const size_t bra_index,
                    const int    ket_angular_momentum,
                    const size_t ket_index) const -> const double *
    {
        return diagonal_values(index) +
               _diagonal_blocks[index].element_offset(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    /// @brief Sets type of matrix.
    /// @param mat_type The type of matrix.
    auto
    set_type(const mat_t mat_type) -> void
    {
        _type = mat_type;
    }

    /// @brief Gets type of matrix.
    /// @return The type of matrix.
    auto
    get_type() const -> mat_t
    {
        return _type;
    }

    /// @brief Gets sparsity patterns of the off-diagonal blocks.
    /// @return The constant reference to vector of sparsity patterns.
    auto
    pair_blocks() const -> const std::vector<CAtomBasisPairSparsity> &
    {
        return _pair_blocks;
    }

    /// @brief Gets sparsity patterns of the diagonal blocks.
    /// @return The constant reference to vector of sparsity patterns.
    auto
    diagonal_blocks() const -> const std::vector<CAtomBasisDiagonalSparsity> &
    {
        return _diagonal_blocks;
    }

    /// @brief Gets sparsity pattern of the off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @return The constant reference to sparsity pattern.
    auto
    pair_block(const size_t index) const -> const CAtomBasisPairSparsity &
    {
        errors::assertMsgCritical(index < _pair_blocks.size(), std::string("SparseMatrix.pair_block: Index of block is out of range"));

        return _pair_blocks[index];
    }

    /// @brief Gets sparsity pattern of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The constant reference to sparsity pattern.
    auto
    diagonal_block(const size_t index) const -> const CAtomBasisDiagonalSparsity &
    {
        errors::assertMsgCritical(index < _diagonal_blocks.size(),
                                  std::string("SparseMatrix.diagonal_block: Index of block is out of range"));

        return _diagonal_blocks[index];
    }

    /// @brief Gets state of the values blocks of matrix.
    /// @return The state of the values blocks.
    auto
    get_values_state() const -> valstat
    {
        return _values_state;
    }

    /// @brief Gets number of values blocks in matrix, i.e. the number of
    /// off-diagonal and diagonal blocks.
    /// @return The number of values blocks.
    auto
    number_of_value_blocks() const -> size_t
    {
        return _values.size();
    }

    /// @brief Gets number of values required to store the integrals of specific
    /// values block.
    /// @param index The index of values block, with the off-diagonal blocks
    /// preceding the diagonal blocks.
    /// @return The number of values.
    auto
    number_of_elements(const size_t index) const -> size_t
    {
        errors::assertMsgCritical(index < _values.size(), std::string("SparseMatrix.number_of_elements: Index of block is out of range"));

        return (index < _pair_blocks.size()) ? _pair_blocks[index].number_of_elements()
                                             : _diagonal_blocks[index - _pair_blocks.size()].number_of_elements();
    }

    /// @brief Gets number of values required to store the integrals of all
    /// values blocks.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        size_t nvals = 0;

        std::ranges::for_each(_pair_blocks, [&](const auto &block) { nvals += block.number_of_elements(); });

        std::ranges::for_each(_diagonal_blocks, [&](const auto &block) { nvals += block.number_of_elements(); });

        return nvals;
    }

    /// @brief Gets memory required to store the values blocks of matrix.
    /// @return The memory in bytes.
    /// @note Only the values blocks are counted, as they are what allocate()
    /// reserves. The sparsity patterns describing the blocks are not included.
    auto
    memory_size() const -> size_t
    {
        return number_of_elements() * sizeof(double);
    }

    /// @brief Gets number of off-diagonal blocks in matrix.
    /// @return The number of off-diagonal blocks.
    auto
    number_of_pair_blocks() const -> size_t
    {
        return _pair_blocks.size();
    }

    /// @brief Gets number of diagonal blocks in matrix.
    /// @return The number of diagonal blocks.
    auto
    number_of_diagonal_blocks() const -> size_t
    {
        return _diagonal_blocks.size();
    }

   private:
    /// @brief Checks that the values blocks are allocated and that a block index
    /// is in range.
    /// @param index The index of block.
    /// @param nblocks The number of blocks of that kind.
    /// @param label The name of the accessor requesting the check.
    auto
    _check_values(const size_t index, const size_t nblocks, const std::string &label) const -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseMatrix.") + label + std::string(": Values blocks of matrix are not allocated"));

        errors::assertMsgCritical(index < nblocks,
                                  std::string("SparseMatrix.") + label + std::string(": Index of block is out of range"));
    }

    /// @brief Deallocates the values blocks of matrix.
    auto
    _deallocate() -> void
    {
        std::ranges::for_each(_values, [](double *block) { delete[] block; });

        std::ranges::fill(_values, nullptr);

        _values_state = valstat::empty;
    }

    /// @brief Allocates the values blocks of matrix and copies the values of
    /// another matrix into them.
    /// @param other The sparse matrix to copy the values of.
    auto
    _copy_values(const CSparseMatrix &other) -> void
    {
        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()), [&](const auto i) {
            const auto nvals = number_of_elements(i);

            _values[i] = new double[nvals];

            std::copy(other._values[i], other._values[i] + nvals, _values[i]);
        });
    }

    /// @brief Adds the sparsity patterns of the non-empty blocks of the atom
    /// basis pair groups.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param groups The atom basis pair groups to describe.
    /// @param screener The integral bound.
    /// @param threshold The screening threshold.
    /// @param storage The storage layout of the diagonal blocks.
    template <typename B>
    auto
    _add_blocks(const CMolecule                  &molecule,
                std::vector<CAtomBasisPairGroup> &groups,
                const B                          &screener,
                const double                      threshold,
                const diagstor                    storage) -> void
    {
        // NOTE: the atom basis pair groups are independent, so their atom pairs are
        // ordered by interatomic distance in parallel. Dynamic scheduling is used
        // as the groups differ widely in the number of atom pairs.

        const auto ngroups = static_cast<int>(groups.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            groups[i].sort_by_distance(molecule);
        }

        // NOTE: the blocks are added in the order of the groups, so that the
        // layout of the values blocks does not depend on the scheduling above.

        std::ranges::for_each(groups, [&](const auto &group) {
            const auto pair_block = CAtomBasisPairSparsity(group, screener, threshold);

            if (pair_block.number_of_pairs() > 0) _pair_blocks.push_back(pair_block);

            const auto diagonal_block = CAtomBasisDiagonalSparsity(group, storage);

            if (diagonal_block.number_of_atoms() > 0) _diagonal_blocks.push_back(diagonal_block);
        });
    }

    /// @brief The sparsity patterns of the off-diagonal blocks.
    std::vector<CAtomBasisPairSparsity> _pair_blocks;

    /// @brief The sparsity patterns of the diagonal blocks.
    std::vector<CAtomBasisDiagonalSparsity> _diagonal_blocks;

    /// @brief The values blocks of matrix, one pointer per block, with the
    /// off-diagonal blocks preceding the diagonal blocks. The pointers are null
    /// while the values blocks are not allocated.
    std::vector<double *> _values;

    /// @brief The type of matrix.
    mat_t _type;

    /// @brief The state of the values blocks of matrix.
    valstat _values_state;
};

#endif /* SparseMatrix_hpp */
