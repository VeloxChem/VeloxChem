//
//  Tabula — custom-recursion molecular-integral machinery.
//  Mixed-precision block-sparse matrix storage.
//

#include "TabulaMixedPrecisionBlockSparseMatrix.hpp"

#include <algorithm>
#include <cmath>

namespace tabula {  // tabula namespace

MixedPrecisionBlockSparseMatrix::MixedPrecisionBlockSparseMatrix()

    : _dimension{0}

    , _groupGlobalAO{}

    , _blocks{}

    , _doubleValues{}

    , _floatValues{}

    , _precisionThreshold{0.0}
{
}

MixedPrecisionBlockSparseMatrix::MixedPrecisionBlockSparseMatrix(const BlockSparseMatrix& source, const double precisionThreshold)

    : _dimension{source.dimension()}

    , _groupGlobalAO{}

    , _blocks{}

    , _doubleValues{}

    , _floatValues{}

    , _precisionThreshold{precisionThreshold}
{
    // carry over the AO group map

    _groupGlobalAO.reserve(source.number_of_groups());

    for (std::size_t g = 0; g < source.number_of_groups(); g++)
    {
        _groupGlobalAO.push_back(source.group_global_ao(g));
    }

    // classify each block by its largest-magnitude element

    _blocks.reserve(source.number_of_blocks());

    for (std::size_t b = 0; b < source.number_of_blocks(); b++)
    {
        const auto src = source.block(b);

        double blockMax = 0.0;

        for (std::size_t i = 0; i < src.rowCount; i++)
        {
            for (std::size_t j = 0; j < src.columnCount; j++)
            {
                blockMax = std::max(blockMax, std::abs(source.value(b, i, j)));
            }
        }

        // a diagonal (same-group) block is always kept double; only an
        // off-diagonal block at or below the threshold demotes to float
        const bool singlePrecision = (src.groupA != src.groupB) && (blockMax <= precisionThreshold);

        const auto offset = singlePrecision ? _floatValues.size() : _doubleValues.size();

        for (std::size_t i = 0; i < src.rowCount; i++)
        {
            for (std::size_t j = 0; j < src.columnCount; j++)
            {
                const auto v = source.value(b, i, j);

                if (singlePrecision)
                {
                    _floatValues.push_back(static_cast<float>(v));
                }
                else
                {
                    _doubleValues.push_back(v);
                }
            }
        }

        _blocks.push_back(Block{src.groupA, src.groupB, src.rowCount, src.columnCount, offset, singlePrecision});
    }
}

auto
MixedPrecisionBlockSparseMatrix::dimension() const -> std::size_t
{
    return _dimension;
}

auto
MixedPrecisionBlockSparseMatrix::number_of_blocks() const -> std::size_t
{
    return _blocks.size();
}

auto
MixedPrecisionBlockSparseMatrix::block(const std::size_t index) const -> Block
{
    return _blocks[index];
}

auto
MixedPrecisionBlockSparseMatrix::group_global_ao(const std::size_t group) const -> std::vector<std::size_t>
{
    return _groupGlobalAO[group];
}

auto
MixedPrecisionBlockSparseMatrix::value(const std::size_t blockIndex, const std::size_t row, const std::size_t column) const -> double
{
    const auto& blk = _blocks[blockIndex];

    const auto index = blk.offset + row * blk.columnCount + column;

    return blk.isSinglePrecision ? static_cast<double>(_floatValues[index]) : _doubleValues[index];
}

auto
MixedPrecisionBlockSparseMatrix::precision_threshold() const -> double
{
    return _precisionThreshold;
}

auto
MixedPrecisionBlockSparseMatrix::single_block_count() const -> std::size_t
{
    std::size_t count = 0;

    for (const auto& blk : _blocks)
    {
        if (blk.isSinglePrecision) count++;
    }

    return count;
}

auto
MixedPrecisionBlockSparseMatrix::double_block_count() const -> std::size_t
{
    return _blocks.size() - single_block_count();
}

auto
MixedPrecisionBlockSparseMatrix::stored_element_count() const -> std::size_t
{
    return _doubleValues.size() + _floatValues.size();
}

auto
MixedPrecisionBlockSparseMatrix::stored_byte_count() const -> std::size_t
{
    return sizeof(double) * _doubleValues.size() + sizeof(float) * _floatValues.size();
}

auto
MixedPrecisionBlockSparseMatrix::to_dense() const -> DenseMatrix
{
    DenseMatrix dense(_dimension, _dimension, Symmetry::symmetric);

    for (std::size_t b = 0; b < _blocks.size(); b++)
    {
        const auto& blk      = _blocks[b];
        const auto& rowAO    = _groupGlobalAO[blk.groupA];
        const auto& columnAO = _groupGlobalAO[blk.groupB];

        for (std::size_t i = 0; i < blk.rowCount; i++)
        {
            for (std::size_t j = 0; j < blk.columnCount; j++)
            {
                const auto v = value(b, i, j);

                dense(rowAO[i], columnAO[j]) = v;

                // an off-diagonal block also fills its transpose
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
