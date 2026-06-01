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
#include <cstdint>
#include <optional>
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

    /// @brief Reserves capacity for the given number of blocks, avoiding rehashes
    /// when the final block count is known up front (e.g. before a bulk merge).
    /// @param count The expected number of blocks.
    auto reserve(const std::size_t count) -> void;

    /// @brief Appends a pre-canonicalized block by copying its payload from a raw
    /// pointer into the data arena. The caller is responsible for canonicalization
    /// (for sym/antisym: i <= j, and a diagonal block already packed lower-
    /// triangular). Used by drivers that build blocks in their own arenas and
    /// merge them in bulk, avoiding a per-block Block allocation.
    /// @param i The (canonical) bra contracted GTO index.
    /// @param j The (canonical) ket contracted GTO index.
    /// @param nrows The number of rows.
    /// @param ncols The number of columns.
    /// @param kind The storage layout (full, or lower_triangular for a packed diagonal).
    /// @param src Pointer to the payload (payload_size(nrows, ncols, kind) values).
    auto add_raw(const int i, const int j, const std::size_t nrows, const std::size_t ncols, const Kind kind, const double *src) -> void;

    /// @brief Sets all block values to zero, keeping the block structure.
    auto zero() -> void;

    /// @brief Gets the symmetry of the matrix.
    /// @return The symmetry of the matrix.
    auto symmetry() const -> SymmetryType;

    /// @brief Checks whether a block exists at the given key.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @return True if a block exists at the key.
    auto contains(const Key &key) const -> bool;

    /// @brief Gets a read-only snapshot of the block at the given key.
    ///
    /// The block data is stored in a contiguous arena, so this returns a copy
    /// (not a reference into storage); mutating it does not affect the matrix.
    /// @param key The (bra, ket) pair of contracted GTO indices.
    /// @return The block, or std::nullopt if absent.
    auto block(const Key &key) const -> std::optional<Block>;

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
    /// @brief Descriptor of one stored block: its canonical (bra, ket) key, its
    /// dimensions and layout, and the offset of its payload in the data arena.
    struct BlockMeta
    {
        int i, j;
        std::uint32_t nrows, ncols;
        std::size_t offset;
        Kind kind;
    };

    /// @brief Ensures _meta is sorted by key with duplicates collapsed (last add
    /// wins). Lazy: add() merely appends, and queries trigger one sort/dedup.
    auto ensure_sorted() const -> void;

    /// @brief The symmetry of the matrix.
    SymmetryType _symmetry;

    /// @brief Per-block descriptors; sorted by key once a query needs ordering.
    mutable std::vector<BlockMeta> _meta;

    /// @brief Whether _meta is currently sorted and deduplicated.
    mutable bool _sorted;

    /// @brief Contiguous arena holding every block's payload back-to-back. One
    /// allocation instead of a heap vector per block.
    std::vector<double> _data;
};

}  // namespace newints

#endif /* newints_SparseMatrix_hpp */
