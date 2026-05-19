//
//  Tabula — custom-recursion molecular-integral machinery.
//  Mixed-precision block-sparse matrix storage.
//

#ifndef TabulaMixedPrecisionBlockSparseMatrix_hpp
#define TabulaMixedPrecisionBlockSparseMatrix_hpp

#include <cstddef>
#include <vector>

#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"

namespace tabula {  // tabula namespace

/// @brief A symmetric block-sparse matrix with mixed precision: each kept
/// group-pair block is held in `double` or `float` according to its
/// magnitude.
///
/// An off-diagonal block whose largest-magnitude element is at or below
/// `precisionThreshold` is stored single-precision, halving that block's
/// footprint; the absolute error it then carries is bounded by
/// `precisionThreshold * 2^-24`. Blocks above the threshold — and every
/// diagonal (same-group) block — stay double-precision.
///
/// This is a storage representation, built from an already-evaluated
/// `BlockSparseMatrix` by classifying each block on its actual magnitude.
class MixedPrecisionBlockSparseMatrix
{
   public:
    /// @brief One stored group-pair block. `offset` indexes the double-value
    /// buffer when `isSinglePrecision` is false, otherwise the float-value
    /// buffer; element `(i,j)` is at `offset + i * columnCount + j`.
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
        /// @brief The block's start in its value buffer.
        std::size_t offset;
        /// @brief Whether the block is stored single-precision.
        bool isSinglePrecision;
    };

    /// @brief Creates an empty mixed-precision block-sparse matrix.
    MixedPrecisionBlockSparseMatrix();

    /// @brief Builds a mixed-precision matrix from a fully-evaluated
    /// block-sparse matrix. Every off-diagonal block whose largest-magnitude
    /// element is at or below `precisionThreshold` is down-converted to
    /// single precision; diagonal (same-group) blocks and above-threshold
    /// blocks stay double.
    /// @param source The fully-evaluated block-sparse matrix.
    /// @param precisionThreshold The block down-conversion threshold.
    MixedPrecisionBlockSparseMatrix(const BlockSparseMatrix& source, const double precisionThreshold);

    /// @brief Gets the dimension of the full matrix.
    /// @return The full-matrix dimension.
    auto dimension() const -> std::size_t;

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

    /// @brief Element access within a block — widened to double.
    /// @param blockIndex The block index.
    /// @param row The block-local row index.
    /// @param column The block-local column index.
    /// @return The value of the addressed element.
    auto value(const std::size_t blockIndex, const std::size_t row, const std::size_t column) const -> double;

    /// @brief Gets the precision threshold the matrix was built with.
    /// @return The precision threshold.
    auto precision_threshold() const -> double;

    /// @brief Gets the number of blocks stored single-precision.
    /// @return The single-precision block count.
    auto single_block_count() const -> std::size_t;

    /// @brief Gets the number of blocks stored double-precision.
    /// @return The double-precision block count.
    auto double_block_count() const -> std::size_t;

    /// @brief Gets the total number of stored scalar elements.
    /// @return The stored element count.
    auto stored_element_count() const -> std::size_t;

    /// @brief Gets the stored footprint in bytes — `8` per double element,
    /// `4` per float element.
    /// @return The stored byte count.
    auto stored_byte_count() const -> std::size_t;

    /// @brief Reconstructs the dense, symmetric matrix in global AO order,
    /// widening the single-precision blocks.
    /// @return The dense matrix.
    auto to_dense() const -> DenseMatrix;

   private:
    /// @brief The dimension of the full matrix.
    std::size_t _dimension;

    /// @brief Per group, its global AO indices in group-local order.
    std::vector<std::vector<std::size_t>> _groupGlobalAO;

    /// @brief The stored blocks.
    std::vector<Block> _blocks;

    /// @brief The concatenated double-precision block data.
    std::vector<double> _doubleValues;

    /// @brief The concatenated single-precision block data.
    std::vector<float> _floatValues;

    /// @brief The precision threshold the matrix was built with.
    double _precisionThreshold;
};

}  // namespace tabula

#endif /* TabulaMixedPrecisionBlockSparseMatrix_hpp */
