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

#ifndef newints_SparseMatrix_hpp
#define newints_SparseMatrix_hpp

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

class CMolecularBasis;
class CDenseMatrix;

namespace newints {

/// @brief Enumerates the symmetry of a sparse matrix.
enum class SymmetryType
{
    symmetric,      ///< matrix equals its transpose
    antisymmetric,  ///< matrix equals the negative of its transpose
    general         ///< no symmetry
};

/// @brief Enumerates the storage layout of a block.
enum class Kind
{
    full,             ///< all nrows * ncols values stored row-major
    lower_triangular  ///< square block; lower triangle (r >= c) packed as n(n+1)/2 values
};

/// @brief A single dense block of a sparse matrix, i.e. the integrals of one
/// pair of contracted GTO shells.
///
/// The block holds the (2 l_a + 1) x (2 l_b + 1) spherical components of the
/// bra (rows) and ket (columns) shells in row-major order.
struct Block
{
    /// @brief The number of rows, i.e. 2 l_a + 1 spherical components of the bra shell.
    std::size_t nrows;

    /// @brief The number of columns, i.e. 2 l_b + 1 spherical components of the ket shell.
    std::size_t ncols;

    /// @brief The block values, row-major. For Kind::full there are nrows * ncols
    /// values; for Kind::lower_triangular the block is square and only its lower
    /// triangle (r >= c) is stored as n(n+1)/2 packed values.
    std::vector<double> data;

    /// @brief The storage layout of the block.
    Kind kind = Kind::full;

    /// @brief The equality operator.
    /// @param other The block to compare with.
    /// @return True if both blocks have identical dimensions, layout and values.
    auto operator==(const Block &other) const -> bool;
};

/// @brief A block-sparse matrix keyed by pairs of contracted GTO indices.
///
/// Keys are (bra contracted GTO index, ket contracted GTO index). The value of
/// each key is a Block of the (2 l_a + 1) x (2 l_b + 1) integrals of that shell
/// pair. The matrix carries a symmetry tag (symmetric, antisymmetric, general).
class SparseMatrix
{
   public:
    /// @brief Key type: the (bra, ket) pair of contracted GTO indices.
    using Key = std::pair<int, int>;

    /// @brief The default constructor, creating a general empty matrix.
    SparseMatrix();

    /// @brief The constructor with a given symmetry.
    /// @param symmetry The symmetry of the matrix.
    explicit SparseMatrix(const SymmetryType symmetry);

    /// @brief The default copy constructor.
    /// @param other The matrix to be copied.
    SparseMatrix(const SparseMatrix &other) = default;

    /// @brief The default move constructor.
    /// @param other The matrix to be moved.
    SparseMatrix(SparseMatrix &&other) noexcept = default;

    /// @brief The default destructor.
    ~SparseMatrix() = default;

    /// @brief The default copy assignment operator.
    /// @param other The matrix to be copy assigned.
    /// @return The assigned matrix.
    auto operator=(const SparseMatrix &other) -> SparseMatrix & = default;

    /// @brief The default move assignment operator.
    /// @param other The matrix to be move assigned.
    /// @return The assigned matrix.
    auto operator=(SparseMatrix &&other) noexcept -> SparseMatrix & = default;

    /// @brief The equality operator.
    /// @param other The matrix to be compared.
    /// @return True if both matrices have identical symmetry and blocks.
    auto operator==(const SparseMatrix &other) const -> bool;

    /// @brief Sets the symmetry of the matrix.
    /// @param symmetry The symmetry to set.
    auto set_symmetry(const SymmetryType symmetry) -> void;

    /// @brief Adds (inserts or overwrites) a block at the given key.
    ///
    /// For a general matrix the full block is stored as given. For a symmetric
    /// or antisymmetric matrix only the upper block triangle (i <= j) is kept:
    /// a block given with i > j is canonicalized to (j, i) by transposing it
    /// (negating it as well when antisymmetric), and a diagonal block (i == j)
    /// is stored packed lower-triangular (Kind::lower_triangular).
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @param block The full block of shell-pair integrals.
    auto add(const Key &key, const Block &block) -> void;

    /// @brief Adds (inserts or overwrites) a block at the given indices.
    /// @param i The bra contracted GTO index.
    /// @param j The ket contracted GTO index.
    /// @param block The full block of shell-pair integrals.
    auto add(const int i, const int j, const Block &block) -> void;

    /// @brief Adds (inserts or overwrites) a block at the given key, moving it.
    ///
    /// Same canonicalization as the copying overload, but for the stored-as-given
    /// cases (general, and i < j for sym/antisym) the block is moved instead of
    /// copied, avoiding a heap allocation of its data vector.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @param block The full block of shell-pair integrals (consumed).
    auto add(const Key &key, Block &&block) -> void;

    /// @brief Adds (inserts or overwrites) a block at the given indices, moving it.
    /// @param i The bra contracted GTO index.
    /// @param j The ket contracted GTO index.
    /// @param block The full block of shell-pair integrals (consumed).
    auto add(const int i, const int j, Block &&block) -> void;

    /// @brief Sets all block values to zero, keeping the block structure.
    auto zero() -> void;

    /// @brief Gets the symmetry of the matrix.
    /// @return The symmetry of the matrix.
    auto symmetry() const -> SymmetryType;

    /// @brief Checks whether a block exists at the given key.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @return True if a block exists at the key.
    auto contains(const Key &key) const -> bool;

    /// @brief Gets a mutable pointer to the block at the given key.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @return The pointer to the block, or nullptr if absent.
    auto block(const Key &key) -> Block *;

    /// @brief Gets a constant pointer to the block at the given key.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @return The constant pointer to the block, or nullptr if absent.
    auto block(const Key &key) const -> const Block *;

    /// @brief Gets the number of stored blocks.
    /// @return The number of stored blocks.
    auto number_of_blocks() const -> std::size_t;

    /// @brief Gets the keys of all stored blocks in ascending order.
    /// @return The vector of keys.
    auto keys() const -> std::vector<Key>;

    /// @brief Converts to a full dense matrix in VeloxChem atomic-orbital
    /// ordering (angular-momentum major: l -> m -> atom -> contraction). For
    /// symmetric / antisymmetric matrices the stored triangle (and packed
    /// diagonal blocks) are expanded into the full matrix.
    /// @param basis The molecular basis defining the indexing and ordering;
    /// must be the basis whose outline produced the stored keys.
    /// @return The dense matrix.
    auto to_dense(const CMolecularBasis &basis) const -> CDenseMatrix;

   private:
    /// @brief The symmetry of the matrix.
    SymmetryType _symmetry;

    /// @brief The stored blocks keyed by (bra, ket) contracted GTO indices.
    std::map<Key, Block> _blocks;
};

}  // namespace newints

#endif /* newints_SparseMatrix_hpp */
