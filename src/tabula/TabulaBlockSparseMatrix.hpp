//
//  Tabula — custom-recursion molecular-integral machinery.
//  Block-sparse matrix storage.
//

#ifndef TabulaBlockSparseMatrix_hpp
#define TabulaBlockSparseMatrix_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "TabulaDenseMatrix.hpp"

namespace tabula {  // tabula namespace

/// @brief A symmetric matrix held in block-sparse storage over an arbitrary
/// atomic-orbital (AO) group partition.
///
/// Only the group-pair blocks handed to the constructor are stored — each a
/// dense `rowCount x columnCount` sub-matrix in group-local AO order, all
/// blocks concatenated in one flat buffer. For an extended molecule the bulk
/// of the matrix is negligible and dropped, giving O(N) storage where a dense
/// matrix is O(N^2).
///
/// A "group" is an opaque partition of the AO space — atoms, shells, whatever
/// the integral driver chooses. `groupGlobalAO[g]` lists group `g`'s global AO
/// indices, mapping a block's local `(i,j)` back to the full matrix.
class BlockSparseMatrix
{
   public:
    /// @brief One stored group-pair block. The block holds the integral of
    /// group `groupA`'s orbitals (rows) with group `groupB`'s (columns),
    /// row-major in group-local AO order; element `(i,j)` is at
    /// `offset + i * columnCount + j` of the value buffer.
    struct Block
    {
        /// @brief The row group index.
        std::size_t groupA;
        /// @brief The column group index.
        std::size_t groupB;
        /// @brief The row count — group A's AO count.
        std::size_t rowCount;
        /// @brief The column count — group B's AO count.
        std::size_t columnCount;
        /// @brief The block's start in the value buffer.
        std::size_t offset;
    };

    /// @brief Creates an empty block-sparse matrix.
    BlockSparseMatrix();

    /// @brief Creates a block-sparse matrix and lays out its (zeroed) storage.
    /// @param dimension The dimension of the full matrix (total AO count).
    /// @param groupGlobalAO Per group, its global AO indices in group-local
    /// order.
    /// @param groupPairs The (groupA, groupB) pairs to store as blocks.
    BlockSparseMatrix(const std::size_t                                       dimension,
                      const std::vector<std::vector<std::size_t>>&            groupGlobalAO,
                      const std::vector<std::pair<std::size_t, std::size_t>>& groupPairs);

    /// @brief Gets the dimension of the full matrix.
    /// @return The full-matrix dimension.
    auto dimension() const -> std::size_t;

    /// @brief Gets the number of AO groups.
    /// @return The number of groups.
    auto number_of_groups() const -> std::size_t;

    /// @brief Gets the number of stored blocks.
    /// @return The number of blocks.
    auto number_of_blocks() const -> std::size_t;

    /// @brief Gets a stored block descriptor.
    /// @param index The block index.
    /// @return The block descriptor.
    auto block(const std::size_t index) const -> Block;

    /// @brief Gets the global AO indices of a group.
    /// @param group The group index.
    /// @return The group's global AO indices.
    auto group_global_ao(const std::size_t group) const -> std::vector<std::size_t>;

    /// @brief Element access within a block.
    /// @param blockIndex The block index.
    /// @param row The block-local row index.
    /// @param column The block-local column index.
    /// @return The value of the addressed element.
    auto value(const std::size_t blockIndex, const std::size_t row, const std::size_t column) const -> double;

    /// @brief Sets an element within a block.
    /// @param blockIndex The block index.
    /// @param row The block-local row index.
    /// @param column The block-local column index.
    /// @param value The value to set.
    auto set_value(const std::size_t blockIndex, const std::size_t row, const std::size_t column, const double value) -> void;

    /// @brief Gets a pointer to the flat, concatenated value buffer.
    /// @return The value buffer pointer.
    auto values() -> double*;

    /// @brief Gets a constant pointer to the flat, concatenated value buffer.
    /// @return The constant value buffer pointer.
    auto values() const -> const double*;

    /// @brief Gets the number of stored scalar elements — the block-sparse
    /// footprint, to compare against the dense `dimension^2`.
    /// @return The stored element count.
    auto stored_element_count() const -> std::size_t;

    /// @brief Reconstructs the dense, symmetric matrix in global AO order.
    /// @return The dense matrix.
    auto to_dense() const -> DenseMatrix;

   private:
    /// @brief The dimension of the full matrix.
    std::size_t _dimension;

    /// @brief Per group, its global AO indices in group-local order.
    std::vector<std::vector<std::size_t>> _groupGlobalAO;

    /// @brief The stored blocks.
    std::vector<Block> _blocks;

    /// @brief The flat, concatenated block value buffer.
    std::vector<double> _values;
};

}  // namespace tabula

#endif /* TabulaBlockSparseMatrix_hpp */
