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


#ifndef SimdVariableMatrix_hpp
#define SimdVariableMatrix_hpp

#include <algorithm>
#include <cstddef>
#include <new>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

/// @brief Class CSimdVariableMatrix stores a two-dimensional array of values with
/// a different number of columns in each row, repeated over a number of blocks,
/// in a form suitable for SIMD operations. Every row is padded, so that it starts
/// at a cache line boundary and can be loaded with aligned SIMD instructions.
/// @note The number of columns of a row is given for the rows of a single block,
/// so the number of rows of a block is the number of given column counts and the
/// number of rows of the matrix is that number times the number of blocks. The
/// rows of every block share the same layout.
class CSimdVariableMatrix
{
   public:
    /// @brief The default constructor.
    CSimdVariableMatrix()

        : _data(nullptr)

        , _columns{}

        , _pitches{}

        , _offsets{0}

        , _blocks(0)
    {
    }

    /// @brief The constructor with number of columns in the rows of a block and
    /// number of blocks.
    /// @param columns The number of columns in each row of a block.
    /// @param blocks The number of blocks in matrix.
    CSimdVariableMatrix(const std::vector<size_t> &columns, const size_t blocks)

        : _data(nullptr)

        , _columns(columns)

        , _pitches{}

        , _offsets{}

        , _blocks(blocks)
    {
        // NOTE: each row is padded on its own, so that every row of every block
        // starts at a cache line boundary. The offsets are the running sums of
        // the padded rows, with the padded size of a block as last element.

        _pitches.reserve(_columns.size());

        _offsets.reserve(_columns.size() + 1);

        _offsets.push_back(0);

        std::ranges::for_each(_columns, [&](const auto ncols) {
            _pitches.push_back(simd::pitch_of(ncols));

            _offsets.push_back(_offsets.back() + _pitches.back());
        });

        _allocate();
    }

    /// @brief The copy constructor.
    /// @param other The matrix to be copied.
    CSimdVariableMatrix(const CSimdVariableMatrix &other)

        : _data(nullptr)

        , _columns(other._columns)

        , _pitches(other._pitches)

        , _offsets(other._offsets)

        , _blocks(other._blocks)
    {
        _allocate();

        if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
    }

    /// @brief The move constructor.
    /// @param other The matrix to be moved.
    CSimdVariableMatrix(CSimdVariableMatrix &&other) noexcept

        : _data(other._data)

        , _columns(std::move(other._columns))

        , _pitches(std::move(other._pitches))

        , _offsets(std::move(other._offsets))

        , _blocks(other._blocks)
    {
        other._data = nullptr;

        other._columns.clear();

        other._pitches.clear();

        other._offsets.assign(1, 0);

        other._blocks = 0;
    }

    /// @brief The destructor.
    ~CSimdVariableMatrix()
    {
        _deallocate();
    }

    /// @brief The copy assignment operator.
    /// @param other The matrix to be copy assigned.
    /// @return The assigned matrix.
    auto
    operator=(const CSimdVariableMatrix &other) -> CSimdVariableMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _columns = other._columns;

            _pitches = other._pitches;

            _offsets = other._offsets;

            _blocks = other._blocks;

            _allocate();

            if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
        }

        return *this;
    }

    /// @brief The move assignment operator.
    /// @param other The matrix to be move assigned.
    /// @return The assigned matrix.
    auto
    operator=(CSimdVariableMatrix &&other) noexcept -> CSimdVariableMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _data = other._data;

            _columns = std::move(other._columns);

            _pitches = std::move(other._pitches);

            _offsets = std::move(other._offsets);

            _blocks = other._blocks;

            other._data = nullptr;

            other._columns.clear();

            other._pitches.clear();

            other._offsets.assign(1, 0);

            other._blocks = 0;
        }

        return *this;
    }

    /// @brief Sets all values of matrix, padding included, to zero.
    auto
    zero() -> void
    {
        if (_data != nullptr) std::fill(_data, _data + number_of_elements(), 0.0);
    }

    /// @brief Gets values of matrix.
    /// @return The pointer to the values of matrix.
    auto
    data() -> double *
    {
        return _data;
    }

    /// @brief Gets values of matrix.
    /// @return The constant pointer to the values of matrix.
    auto
    data() const -> const double *
    {
        return _data;
    }

    /// @brief Gets values of specific row of specific block of matrix.
    /// @param block The index of block.
    /// @param row The index of row in block.
    /// @return The pointer to the values of row, aligned to a cache line boundary.
    auto
    data(const size_t block, const size_t row) -> double *
    {
        _check_index(block, row);

        return _data + block * block_pitch() + _offsets[row];
    }

    /// @brief Gets values of specific row of specific block of matrix.
    /// @param block The index of block.
    /// @param row The index of row in block.
    /// @return The constant pointer to the values of row, aligned to a cache line
    /// boundary.
    auto
    data(const size_t block, const size_t row) const -> const double *
    {
        _check_index(block, row);

        return _data + block * block_pitch() + _offsets[row];
    }

    /// @brief Gets number of columns in the rows of a block.
    /// @return The constant reference to vector of number of columns.
    auto
    columns() const -> const std::vector<size_t> &
    {
        return _columns;
    }

    /// @brief Gets padded number of columns in the rows of a block.
    /// @return The constant reference to vector of padded number of columns.
    auto
    pitches() const -> const std::vector<size_t> &
    {
        return _pitches;
    }

    /// @brief Gets number of columns in specific row of a block.
    /// @param row The index of row in block.
    /// @return The number of columns.
    auto
    number_of_columns(const size_t row) const -> size_t
    {
        if (row >= _columns.size())
        {
            errors::assertMsgCritical(false, std::string("SimdVariableMatrix.number_of_columns: Index of row is out of range"));
        }

        return _columns[row];
    }

    /// @brief Gets padded number of columns in specific row of a block.
    /// @param row The index of row in block.
    /// @return The padded number of columns.
    auto
    pitch(const size_t row) const -> size_t
    {
        if (row >= _pitches.size())
        {
            errors::assertMsgCritical(false, std::string("SimdVariableMatrix.pitch: Index of row is out of range"));
        }

        return _pitches[row];
    }

    /// @brief Gets padded number of values in a block of matrix.
    /// @return The padded number of values in a block.
    auto
    block_pitch() const -> size_t
    {
        return _offsets.back();
    }

    /// @brief Gets number of blocks in matrix.
    /// @return The number of blocks.
    auto
    number_of_blocks() const -> size_t
    {
        return _blocks;
    }

    /// @brief Gets number of rows in a block of matrix.
    /// @return The number of rows in a block.
    auto
    number_of_rows_in_block() const -> size_t
    {
        return _columns.size();
    }

    /// @brief Gets number of rows in matrix.
    /// @return The number of rows.
    auto
    number_of_rows() const -> size_t
    {
        return _blocks * _columns.size();
    }

    /// @brief Gets number of values in matrix, padding included.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        return _blocks * block_pitch();
    }

    /// @brief Gets memory required to store the values of matrix.
    /// @return The memory in bytes.
    auto
    memory_size() const -> size_t
    {
        return number_of_elements() * sizeof(double);
    }

   private:
    /// @brief Checks that a block and row index are in range.
    /// @param block The index of block.
    /// @param row The index of row in block.
    auto
    _check_index(const size_t block, const size_t row) const -> void
    {
        // NOTE: the messages are built only when a check fails, as constructing
        // them eagerly allocates on every access of a row.

        if (block >= _blocks)
        {
            errors::assertMsgCritical(false, std::string("SimdVariableMatrix.data: Index of block is out of range"));
        }

        if (row >= _columns.size())
        {
            errors::assertMsgCritical(false, std::string("SimdVariableMatrix.data: Index of row is out of range"));
        }
    }

    /// @brief Allocates the values of matrix, leaving their content undefined.
    auto
    _allocate() -> void
    {
        if (const auto nelems = number_of_elements(); nelems > 0)
        {
            _data = static_cast<double *>(::operator new[](nelems * sizeof(double), std::align_val_t{simd::cache_line_size()}));
        }
    }

    /// @brief Deallocates the values of matrix.
    auto
    _deallocate() -> void
    {
        if (_data != nullptr)
        {
            ::operator delete[](_data, std::align_val_t{simd::cache_line_size()});

            _data = nullptr;
        }
    }

    /// @brief The values of matrix, stored block wise with padded rows.
    double *_data;

    /// @brief The number of columns in each row of a block.
    std::vector<size_t> _columns;

    /// @brief The padded number of columns in each row of a block.
    std::vector<size_t> _pitches;

    /// @brief The offsets of the rows of a block, with the padded number of
    /// values in a block as last element.
    std::vector<size_t> _offsets;

    /// @brief The number of blocks in matrix.
    size_t _blocks;
};

#endif /* SimdVariableMatrix_hpp */
