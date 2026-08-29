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
#include <vector>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"

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
class CSparseMatrix
{
   public:
    /// @brief The default constructor.
    CSparseMatrix()

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _type(mat_t::general)
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

        , _type(mat_type)
    {
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
    /// @brief The sparsity patterns of the off-diagonal blocks.
    std::vector<CAtomBasisPairSparsity> _pair_blocks;

    /// @brief The sparsity patterns of the diagonal blocks.
    std::vector<CAtomBasisDiagonalSparsity> _diagonal_blocks;

    /// @brief The type of matrix.
    mat_t _type;
};

#endif /* SparseMatrix_hpp */
