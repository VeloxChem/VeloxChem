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

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "DenseMatrix.hpp"
#include "Matrix.hpp"
#include "ScreenedBasisFunctionPair.hpp"

/// @brief Class CSparseMatrix stores a block sparse matrix whose block layout is
/// defined by a vector of screened basis function pairs. Each block is a dense
/// matrix of (2 l_a + 1) (2 l_b + 1) rows and npairs columns, where l_a and l_b
/// are the angular momenta of the bra and ket basis functions of the screened
/// pair and npairs is the number of surviving atom pairs. Blocks are keyed by
/// the index of the screened basis function pair and are created lazily, sized
/// from the screened pairs supplied at construction. The matrix carries a
/// symmetry type (general, symmetric or antisymmetric).
class CSparseMatrix
{
   public:
    /// @brief The default constructor.
    CSparseMatrix();

    /// @brief The constructor sizing the block layout from screened basis
    /// function pairs.
    /// @param screened_pairs The vector of screened basis function pairs.
    /// @param mtype The symmetry type of the matrix.
    CSparseMatrix(const std::vector<CScreenedBasisFunctionPair> &screened_pairs, const mat_t mtype);

    /// @brief The default copy constructor.
    /// @param other The sparse matrix to be copied.
    CSparseMatrix(const CSparseMatrix &other);

    /// @brief The default move constructor.
    /// @param other The sparse matrix to be moved.
    CSparseMatrix(CSparseMatrix &&other) noexcept;

    /// @brief The default destructor.
    ~CSparseMatrix() = default;

    /// @brief The default copy assignment operator.
    /// @param other The sparse matrix to be copy assigned.
    /// @return The assigned sparse matrix.
    auto operator=(const CSparseMatrix &other) -> CSparseMatrix &;

    /// @brief The default move assignment operator.
    /// @param other The sparse matrix to be move assigned.
    /// @return The assigned sparse matrix.
    auto operator=(CSparseMatrix &&other) noexcept -> CSparseMatrix &;

    /// @brief The equality operator.
    /// @param other The sparse matrix to be compared.
    /// @return True if sparse matrices are equal, False otherwise.
    auto operator==(const CSparseMatrix &other) const -> bool;

    /// @brief The non-equality operator.
    /// @param other The sparse matrix to be compared.
    /// @return True if sparse matrices are not equal, False otherwise.
    auto operator!=(const CSparseMatrix &other) const -> bool;

    /// @brief Gets the symmetry type of the matrix.
    /// @return The symmetry type.
    auto type() const -> mat_t;

    /// @brief Gets the number of blocks in the block layout (one per screened
    /// basis function pair), whether or not they have been created.
    /// @return The number of blocks in the layout.
    auto number_of_keys() const -> size_t;

    /// @brief Gets the number of currently allocated blocks.
    /// @return The number of allocated blocks.
    auto number_of_blocks() const -> size_t;

    /// @brief Checks whether a block has been allocated for a key.
    /// @param key The index of the screened basis function pair.
    /// @return True if the block is allocated, False otherwise.
    auto has_block(const size_t key) const -> bool;

    /// @brief Gets the number of rows of the block for a key.
    /// @param key The index of the screened basis function pair.
    /// @return The number of rows, (2 l_a + 1) (2 l_b + 1).
    auto block_rows(const size_t key) const -> int;

    /// @brief Gets the number of columns of the block for a key.
    /// @param key The index of the screened basis function pair.
    /// @return The number of columns, npairs.
    auto block_columns(const size_t key) const -> int;

    /// @brief Gets the block for a key, creating a zero block of the layout size
    /// if it has not been allocated yet.
    /// @param key The index of the screened basis function pair.
    /// @return The block of values.
    auto block(const size_t key) -> CDenseMatrix &;

    /// @brief Gets the allocated block for a key.
    /// @param key The index of the screened basis function pair.
    /// @return The block of values.
    auto block(const size_t key) const -> const CDenseMatrix &;

    /// @brief Gets the sorted keys of the allocated blocks.
    /// @return The vector of keys.
    auto keys() const -> std::vector<size_t>;

    /// @brief Zeroes the values of all allocated blocks.
    auto zero() -> void;

   private:
    /// @brief The symmetry type of the matrix.
    mat_t _mtype;

    /// @brief The number of rows of each block in the layout, indexed by key.
    std::vector<int> _block_rows;

    /// @brief The number of columns of each block in the layout, indexed by key.
    std::vector<int> _block_cols;

    /// @brief The allocated blocks of values, keyed by screened basis function
    /// pair index.
    std::unordered_map<size_t, CDenseMatrix> _blocks;
};

#endif /* SparseMatrix_hpp */
