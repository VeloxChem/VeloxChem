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
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairGroup.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "OpenMPFunc.hpp"
#include "TensorComponents.hpp"
#include "ValuesState.hpp"

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
    /// @brief The number of blocks per thread aimed at when the target number of
    /// atom pairs of a block is chosen. The blocks are a few per thread, so that
    /// dynamic scheduling has enough of them to even out the ones which differ in
    /// cost, and no more, as a block carries a fixed cost and the blocks contend
    /// for the memory. Measured on fourteen threads, where two per thread is
    /// five percent better than four and four is twice as good as sixteen.
    static constexpr size_t blocks_per_thread = 2;

    /// @brief The smallest target number of atom pairs of a block chosen. A
    /// block carries a fixed cost which does not shrink with the atom pairs it
    /// holds, chiefly the bisection of the screening over the pairs of
    /// primitives, so a molecule too small to fill the threads is divided into
    /// fewer blocks rather than into blocks whose fixed cost outweighs their work.
    static constexpr size_t min_block_size = 2048;

    /// @brief The number of elements of the dense matrix one thread sets to zero
    /// at a time when the matrix is reconstructed.
    static constexpr size_t _dense_chunk_size = 1 << 20;

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
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    template <typename B>
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage,
                  const size_t           block_size = 0)

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

        _add_blocks(molecule, groups, screener, threshold, storage, block_size);
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
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    template <typename B>
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &bra_basis,
                  const CMolecularBasis &ket_basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage,
                  const size_t           block_size = 0)

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        errors::assertMsgCritical(mat_type == mat_t::general,
                                  std::string("SparseMatrix: Matrix with two molecular bases must be of general type"));

        auto groups = bra_basis.basis_pair_groups(ket_basis);

        _add_blocks(molecule, groups, screener, threshold, storage, block_size);
    }

    /// @brief The constructor with molecule, molecular basis, named integral
    /// bound and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param bound The integral bound to screen atom pairs with.
    /// @param threshold The screening threshold.
    /// @param mat_type The type of matrix.
    /// @param storage The storage layout of the diagonal blocks.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &basis,
                  const screener         bound,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage,
                  const size_t           block_size = 0)

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        auto groups = (mat_type == mat_t::general) ? basis.basis_pair_groups(basis) : basis.basis_pair_groups();

        _add_named_blocks(molecule, groups, bound, threshold, storage, block_size);
    }

    /// @brief The constructor with molecule, molecular bases on bra and ket
    /// sides, named integral bound and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param bound The integral bound to screen atom pairs with.
    /// @param threshold The screening threshold.
    /// @param mat_type The type of matrix, which must be general.
    /// @param storage The storage layout of the diagonal blocks.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    CSparseMatrix(const CMolecule       &molecule,
                  const CMolecularBasis &bra_basis,
                  const CMolecularBasis &ket_basis,
                  const screener         bound,
                  const double           threshold,
                  const mat_t            mat_type,
                  const diagstor         storage,
                  const size_t           block_size = 0)

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        errors::assertMsgCritical(mat_type == mat_t::general,
                                  std::string("SparseMatrix: Matrix with two molecular bases must be of general type"));

        auto groups = bra_basis.basis_pair_groups(ket_basis);

        _add_named_blocks(molecule, groups, bound, threshold, storage, block_size);
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
    /// @note A matrix reaches hundreds of gigabytes, so an allocation which
    /// throws is an outcome to handle rather than one to ignore. The blocks
    /// already allocated are freed before the exception leaves, so the matrix is
    /// left empty and consistent: without that the pointers would survive behind
    /// a state of empty, and the next call would overwrite and leak them.
    auto
    allocate() -> void
    {
        if (_values_state == valstat::allocated) return;

        try
        {
            std::ranges::for_each(std::views::iota(size_t{0}, _values.size()),
                                  [&](const auto i) { _values[i] = new double[number_of_elements(i)]; });
        }
        catch (...)
        {
            _deallocate();

            throw;
        }

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
        _check_values(index, _pair_blocks.size(), "pair_values");

        return _values[index];
    }

    /// @brief Gets values of the off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @return The constant pointer to the values of the block.
    auto
    pair_values(const size_t index) const -> const double *
    {
        _check_values(index, _pair_blocks.size(), "pair_values");

        return _values[index];
    }

    /// @brief Gets values of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The pointer to the values of the block.
    auto
    diagonal_values(const size_t index) -> double *
    {
        _check_values(index, _diagonal_blocks.size(), "diagonal_values");

        return _values[_pair_blocks.size() + index];
    }

    /// @brief Gets values of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The constant pointer to the values of the block.
    auto
    diagonal_values(const size_t index) const -> const double *
    {
        _check_values(index, _diagonal_blocks.size(), "diagonal_values");

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

    /// @brief Reconstructs the dense matrix in atomic orbital basis from the
    /// values blocks of matrix.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param values The values of the dense matrix, as a row major array of
    /// bra_basis.dimensions_of_basis() rows and ket_basis.dimensions_of_basis()
    /// columns, which must be allocated by the caller.
    /// @note The atomic orbitals are ordered as in CMatrix, i.e. by ascending
    /// angular momentum, and within an angular momentum by angular component and
    /// then by basis function. The elements not covered by the sparsity patterns
    /// are set to zero.
    auto
    to_dense(const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis, double *values) const -> void
    {
        errors::assertMsgCritical(_values_state != valstat::empty,
                                  std::string("SparseMatrix.to_dense: Values blocks of matrix are not allocated"));

        const auto nrows = bra_basis.dimensions_of_basis();

        const auto ncols = ket_basis.dimensions_of_basis();

        // NOTE: the elements are set to zero by the threads, as the dense matrix
        // reaches several gigabytes and its zeroing is a large part of the
        // reconstruction. Static scheduling is used, as the chunks cost the same
        // and each thread then zeroes a contiguous stretch of the matrix.

        const auto nchunks = static_cast<int>((nrows * ncols + _dense_chunk_size - 1) / _dense_chunk_size);

#pragma omp parallel for schedule(static) if (nchunks > 1)
        for (int i = 0; i < nchunks; i++)
        {
            const auto first = static_cast<size_t>(i) * _dense_chunk_size;

            const auto last = std::min(first + _dense_chunk_size, nrows * ncols);

            std::fill(values + first, values + last, 0.0);
        }

        // NOTE: the dense index of an atomic orbital is looked up instead of
        // recomputed, as the innermost loops below run over the atom pairs and
        // would otherwise scan the basis for every value.

        const auto bra_starts = denseidx::make_dense_starts(bra_basis);

        const auto ket_starts = denseidx::make_dense_starts(ket_basis);

        const auto bra_strides = denseidx::make_dense_strides(bra_basis);

        const auto ket_strides = denseidx::make_dense_strides(ket_basis);

        const auto bra_nmoms = bra_strides.size();

        const auto ket_nmoms = ket_strides.size();

        // NOTE: an element of the dense matrix belongs to a single atom pair, and
        // an atom pair belongs to a single block, so the blocks write to disjoint
        // elements and need no synchronization. Dynamic scheduling is used as the
        // blocks differ in the number of the combinations of basis functions.

        const auto nblocks = static_cast<int>(_pair_blocks.size());

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            const auto &block = _pair_blocks[static_cast<size_t>(iblk)];

            const auto *block_values = pair_values(static_cast<size_t>(iblk));

            const auto &bra_atoms = block.bra_atoms();

            const auto &ket_atoms = block.ket_atoms();

            const auto a_indices = denseidx::index_functions(bra_basis.basis_set(block.bra_index()));

            const auto b_indices = denseidx::index_functions(ket_basis.basis_set(block.ket_index()));

            for (size_t i = 0; i < a_indices.size(); i++)
            {
                const auto [la, ia] = a_indices[i];

                for (size_t j = 0; j < b_indices.size(); j++)
                {
                    const auto [lb, jb] = b_indices[j];

                    const auto npairs = block.number_of_pairs(la, ia, lb, jb);

                    if (npairs == 0) continue;

                    const auto *cell = block_values + block.element_offset(la, ia, lb, jb);

                    const auto bra_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

                    const auto ket_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

                    // NOTE: the values of a combination of basis functions are
                    // stored as one row of atom pairs per angular component, with
                    // the component on bra side as the slowest running index.

                    for (size_t ma = 0; ma < bra_ncomps; ma++)
                    {
                        for (size_t mb = 0; mb < ket_ncomps; mb++)
                        {
                            const auto *prim = cell + (ma * ket_ncomps + mb) * npairs;

                            for (size_t k = 0; k < npairs; k++)
                            {
                                const auto row = bra_starts[static_cast<size_t>(bra_atoms[k]) * bra_nmoms + la] + ia + ma * bra_strides[la];

                                const auto col = ket_starts[static_cast<size_t>(ket_atoms[k]) * ket_nmoms + lb] + jb + mb * ket_strides[lb];

                                values[row * ncols + col] = prim[k];

                                // NOTE: only one of the two atom pairs is stored
                                // when the matrix is symmetric or antisymmetric,
                                // so the transposed element is filled here.

                                if (_type == mat_t::symmetric) values[col * ncols + row] = prim[k];

                                if (_type == mat_t::antisymmetric) values[col * ncols + row] = -prim[k];
                            }
                        }
                    }
                }
            }
        }

        _diagonal_to_dense(bra_basis, ket_basis, bra_starts, ket_starts, bra_strides, ket_strides, values, ncols);
    }

    /// @brief Reconstructs the dense matrix in atomic orbital basis from the
    /// values blocks of matrix.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param values The values of the dense matrix, as a row major array of
    /// basis.dimensions_of_basis() rows and columns, which must be allocated by
    /// the caller.
    auto
    to_dense(const CMolecularBasis &basis, double *values) const -> void
    {
        to_dense(basis, basis, values);
    }

   private:
    /// @brief Adds the values of the diagonal blocks to the dense matrix.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param bra_starts The dense indices of the basis functions on bra side.
    /// @param ket_starts The dense indices of the basis functions on ket side.
    /// @param bra_strides The distances between angular components on bra side.
    /// @param ket_strides The distances between angular components on ket side.
    /// @param values The values of the dense matrix.
    /// @param ncols The number of columns of the dense matrix.
    auto
    _diagonal_to_dense(const CMolecularBasis     &bra_basis,
                       const CMolecularBasis     &ket_basis,
                       const std::vector<size_t> &bra_starts,
                       const std::vector<size_t> &ket_starts,
                       const std::vector<size_t> &bra_strides,
                       const std::vector<size_t> &ket_strides,
                       double                    *values,
                       const size_t               ncols) const -> void
    {
        const auto bra_nmoms = bra_strides.size();

        const auto ket_nmoms = ket_strides.size();

        // NOTE: the diagonal blocks carry atoms of their own and write to
        // disjoint elements, as the off-diagonal blocks above do.

        const auto nblocks = static_cast<int>(_diagonal_blocks.size());

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            const auto &block = _diagonal_blocks[static_cast<size_t>(iblk)];

            errors::assertMsgCritical(block.get_storage() == diagstor::scalar,
                                      std::string("SparseMatrix.to_dense: Storage layout of the diagonal blocks is not supported"));

            const auto *block_values = diagonal_values(static_cast<size_t>(iblk));

            const auto &atoms = block.atoms();

            const auto a_indices = denseidx::index_functions(bra_basis.basis_set(block.bra_index()));

            const auto b_indices = denseidx::index_functions(ket_basis.basis_set(block.ket_index()));

            for (size_t i = 0; i < a_indices.size(); i++)
            {
                const auto [la, ia] = a_indices[i];

                for (size_t j = 0; j < b_indices.size(); j++)
                {
                    // NOTE: only the stored combinations are visited, as the
                    // values of the reverse order share their storage.

                    if (block.is_triangular() && (i > j)) continue;

                    const auto [lb, jb] = b_indices[j];

                    if (block.number_of_elements(la, ia, lb, jb) == 0) continue;

                    // NOTE: the operator is spherically symmetric about the atom,
                    // so a single value is stored for the whole block and it is
                    // diagonal in the angular components.

                    const auto fval = block_values[block.element_offset(la, ia, lb, jb)];

                    const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

                    for (size_t k = 0; k < atoms.size(); k++)
                    {
                        const auto atom = static_cast<size_t>(atoms[k]);

                        for (size_t m = 0; m < ncomps; m++)
                        {
                            const auto row = bra_starts[atom * bra_nmoms + la] + ia + m * bra_strides[la];

                            const auto col = ket_starts[atom * ket_nmoms + lb] + jb + m * ket_strides[lb];

                            values[row * ncols + col] = fval;

                            if (block.is_triangular() && (row != col)) values[col * ncols + row] = fval;
                        }
                    }
                }
            }
        }
    }

    /// @brief Adds the blocks screened with a named two-center integral bound.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param groups The atom basis pair groups to describe.
    /// @param bound The integral bound to screen atom pairs with.
    /// @param threshold The screening threshold.
    /// @param storage The storage layout of the diagonal blocks.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    auto
    _add_named_blocks(const CMolecule                  &molecule,
                      std::vector<CAtomBasisPairGroup> &groups,
                      const screener                    bound,
                      const double                      threshold,
                      const diagstor                    storage,
                      const size_t                      block_size) -> void
    {
        if (bound == screener::overlap)
        {
            _add_blocks(molecule, groups, screenfunc::two_center_overlap_bound, threshold, storage, block_size);
        }
        else if (bound == screener::kinetic_energy)
        {
            _add_blocks(molecule, groups, screenfunc::two_center_kinetic_energy_bound, threshold, storage, block_size);
        }
        else if (bound == screener::nuclear_potential)
        {
            _add_blocks(molecule, groups, screenfunc::two_center_nuclear_potential_bound, threshold, storage, block_size);
        }
        else
        {
            errors::assertMsgCritical(false, std::string("SparseMatrix: Integral bound is not a two-center bound"));
        }
    }

    /// @brief Checks that the values blocks are allocated and that a block index
    /// is in range.
    /// @param index The index of block.
    /// @param nblocks The number of blocks of that kind.
    /// @param label The name of the accessor requesting the check.
    auto
    _check_values(const size_t index, const size_t nblocks, const char *label) const -> void
    {
        // NOTE: the message of a check is formed only where the check fails. The
        // accessors guarded here are called once for every combination of basis
        // functions of every block, and forming three messages on every call,
        // which the arguments of the assertions would do, costs more than the
        // checks themselves.

        if (_values_state != valstat::allocated)
        {
            errors::assertMsgCritical(
                false, std::string("SparseMatrix.") + label + std::string(": Values blocks of matrix are not allocated"));
        }

        if (index >= nblocks)
        {
            errors::assertMsgCritical(false, std::string("SparseMatrix.") + label + std::string(": Index of block is out of range"));
        }

        // NOTE: the values blocks are indexed with the off-diagonal blocks first,
        // so the index is checked against them as well as against the blocks of
        // its own kind.

        if (_values.size() != _pair_blocks.size() + _diagonal_blocks.size())
        {
            errors::assertMsgCritical(
                false, std::string("SparseMatrix.") + label + std::string(": Values blocks are inconsistent with blocks"));
        }
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
    /// @note The blocks already allocated are freed before an exception leaves,
    /// as this runs from the copy constructor as well, where a matrix whose
    /// construction did not complete is never destroyed and would leak all of
    /// them.
    auto
    _copy_values(const CSparseMatrix &other) -> void
    {
        try
        {
            std::ranges::for_each(std::views::iota(size_t{0}, _values.size()), [&](const auto i) {
                const auto nvals = number_of_elements(i);

                _values[i] = new double[nvals];

                std::copy(other._values[i], other._values[i] + nvals, _values[i]);
            });
        }
        catch (...)
        {
            _deallocate();

            throw;
        }
    }

    /// @brief Adds the sparsity patterns of the non-empty blocks of the atom
    /// basis pair groups.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param groups The atom basis pair groups to describe.
    /// @param screener The integral bound.
    /// @param threshold The screening threshold.
    /// @param storage The storage layout of the diagonal blocks.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    template <typename B>
    auto
    _add_blocks(const CMolecule                  &molecule,
                std::vector<CAtomBasisPairGroup> &groups,
                const B                          &screener,
                const double                      threshold,
                const diagstor                    storage,
                const size_t                      block_size) -> void
    {
        // NOTE: the atom basis pair groups are as many as the pairs of the
        // unique atom bases, so their number is set by the variety of the
        // elements of the molecule and not by its size, and the largest of them
        // holds a third of the atom pairs. Dividing them into blocks of a target
        // number of atom pairs makes the number of the blocks follow the size of
        // the molecule instead, so the work of every stage below divides for any
        // number of threads.

        const auto nblock_pairs = (block_size == 0) ? CAtomBasisPairGroup::make_block_size(groups, blocks_per_thread, min_block_size) : block_size;

        auto blocks = (nblock_pairs == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, nblock_pairs);

        // NOTE: the atom pairs of all the blocks are ordered by interatomic
        // distance before the sparsity patterns are described, as the patterns
        // are read off the leading atom pairs which survive the screening. A
        // block holds a subrange of the atom pairs of its group and is ordered
        // within itself, which is all the bisection of the screening needs, as
        // the screening keeps an atom pair or drops it on its own distance.

        CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

        // NOTE: the patterns are held in a vector indexed by the block, so that
        // they are described in any order and added in the order of the blocks.
        // The layout of the values blocks therefore does not depend on the
        // scheduling.

        const auto nblocks = static_cast<int>(blocks.size());

        std::vector<std::optional<CAtomBasisPairSparsity>> pair_blocks(blocks.size());

        std::vector<std::optional<CAtomBasisDiagonalSparsity>> diagonal_blocks(blocks.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nblocks; i++)
        {
            pair_blocks[i].emplace(blocks[i], screener, threshold);

            diagonal_blocks[i].emplace(blocks[i], storage);
        }

        // NOTE: an empty pattern is left out, as it carries no values block of
        // its own.

        for (int i = 0; i < nblocks; i++)
        {
            if (pair_blocks[i]->number_of_pairs() > 0) _pair_blocks.push_back(std::move(*pair_blocks[i]));

            if (diagonal_blocks[i]->number_of_atoms() > 0) _diagonal_blocks.push_back(std::move(*diagonal_blocks[i]));
        }

        // NOTE: one values block per sparsity pattern, sized here so that every
        // construction path leaves the values blocks consistent with the blocks.

        _values.assign(_pair_blocks.size() + _diagonal_blocks.size(), nullptr);
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
