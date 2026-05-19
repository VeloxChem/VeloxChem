//
//  Tabula — custom-recursion molecular-integral machinery.
//  Block-sparse matrix storage.
//

#include "TabulaBlockSparseMatrix.hpp"

#include <string>

#include "ErrorHandler.hpp"

namespace tabula {  // tabula namespace

BlockSparseMatrix::BlockSparseMatrix()

    : _dimension{0}

    , _groupGlobalAO{}

    , _blocks{}

    , _values{}
{
}

BlockSparseMatrix::BlockSparseMatrix(const std::size_t                                       dimension,
                                     const std::vector<std::vector<std::size_t>>&            groupGlobalAO,
                                     const std::vector<std::pair<std::size_t, std::size_t>>& groupPairs)

    : _dimension{dimension}

    , _groupGlobalAO{groupGlobalAO}

    , _blocks{}

    , _values{}
{
    _blocks.reserve(groupPairs.size());

    std::size_t offset = 0;

    for (const auto& [groupA, groupB] : groupPairs)
    {
        errors::assertMsgCritical(groupA < groupGlobalAO.size() && groupB < groupGlobalAO.size(),
                                  std::string("BlockSparseMatrix: group index out of range"));

        const auto rowCount    = groupGlobalAO[groupA].size();
        const auto columnCount = groupGlobalAO[groupB].size();

        _blocks.push_back(Block{groupA, groupB, rowCount, columnCount, offset});

        offset += rowCount * columnCount;
    }

    _values = std::vector<double>(offset, 0.0);
}

auto
BlockSparseMatrix::dimension() const -> std::size_t
{
    return _dimension;
}

auto
BlockSparseMatrix::number_of_groups() const -> std::size_t
{
    return _groupGlobalAO.size();
}

auto
BlockSparseMatrix::number_of_blocks() const -> std::size_t
{
    return _blocks.size();
}

auto
BlockSparseMatrix::block(const std::size_t index) const -> Block
{
    return _blocks[index];
}

auto
BlockSparseMatrix::group_global_ao(const std::size_t group) const -> std::vector<std::size_t>
{
    return _groupGlobalAO[group];
}

auto
BlockSparseMatrix::value(const std::size_t blockIndex, const std::size_t row, const std::size_t column) const -> double
{
    const auto& blk = _blocks[blockIndex];

    return _values[blk.offset + row * blk.columnCount + column];
}

auto
BlockSparseMatrix::set_value(const std::size_t blockIndex, const std::size_t row, const std::size_t column, const double value) -> void
{
    const auto& blk = _blocks[blockIndex];

    _values[blk.offset + row * blk.columnCount + column] = value;
}

auto
BlockSparseMatrix::values() -> double*
{
    return _values.data();
}

auto
BlockSparseMatrix::values() const -> const double*
{
    return _values.data();
}

auto
BlockSparseMatrix::stored_element_count() const -> std::size_t
{
    return _values.size();
}

auto
BlockSparseMatrix::to_dense() const -> DenseMatrix
{
    DenseMatrix dense(_dimension, _dimension, Symmetry::symmetric);

    for (const auto& blk : _blocks)
    {
        const auto& rowAO    = _groupGlobalAO[blk.groupA];
        const auto& columnAO = _groupGlobalAO[blk.groupB];

        for (std::size_t i = 0; i < blk.rowCount; i++)
        {
            for (std::size_t j = 0; j < blk.columnCount; j++)
            {
                const auto v = _values[blk.offset + i * blk.columnCount + j];

                dense(rowAO[i], columnAO[j]) = v;

                // an off-diagonal block also fills its transpose, since the
                // stored block set covers only one side of the symmetry
                if (blk.groupA != blk.groupB)
                {
                    dense(columnAO[j], rowAO[i]) = v;
                }
            }
        }
    }

    return dense;
}

}  // namespace tabula
