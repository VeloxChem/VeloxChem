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


#ifndef SimdMatrix_hpp
#define SimdMatrix_hpp

#include <algorithm>
#include <cstddef>
#include <new>
#include <string>
#include <utility>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

/// @brief Class CSimdMatrix stores a two-dimensional array of values in a form
/// suitable for SIMD operations. The rows are padded, so that every row starts
/// at a cache line boundary and can be loaded with aligned SIMD instructions.
class CSimdMatrix
{
   public:
    /// @brief The default constructor.
    CSimdMatrix()

        : _data(nullptr)

        , _rows(0)

        , _columns(0)

        , _pitch(0)
    {
    }

    /// @brief The constructor with number of rows and columns.
    /// @param rows The number of rows in matrix.
    /// @param columns The number of columns in matrix.
    CSimdMatrix(const size_t rows, const size_t columns)

        : _data(nullptr)

        , _rows(rows)

        , _columns(columns)

        , _pitch(simd::pitch_of(columns))
    {
        _allocate();
    }

    /// @brief The copy constructor.
    /// @param other The matrix to be copied.
    CSimdMatrix(const CSimdMatrix &other)

        : _data(nullptr)

        , _rows(other._rows)

        , _columns(other._columns)

        , _pitch(other._pitch)
    {
        _allocate();

        if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
    }

    /// @brief The move constructor.
    /// @param other The matrix to be moved.
    CSimdMatrix(CSimdMatrix &&other) noexcept

        : _data(other._data)

        , _rows(other._rows)

        , _columns(other._columns)

        , _pitch(other._pitch)
    {
        other._data = nullptr;

        other._rows = 0;

        other._columns = 0;

        other._pitch = 0;
    }

    /// @brief The destructor.
    ~CSimdMatrix()
    {
        _deallocate();
    }

    /// @brief The copy assignment operator.
    /// @param other The matrix to be copy assigned.
    /// @return The assigned matrix.
    auto
    operator=(const CSimdMatrix &other) -> CSimdMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _rows = other._rows;

            _columns = other._columns;

            _pitch = other._pitch;

            _allocate();

            if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
        }

        return *this;
    }

    /// @brief The move assignment operator.
    /// @param other The matrix to be move assigned.
    /// @return The assigned matrix.
    auto
    operator=(CSimdMatrix &&other) noexcept -> CSimdMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _data = other._data;

            _rows = other._rows;

            _columns = other._columns;

            _pitch = other._pitch;

            other._data = nullptr;

            other._rows = 0;

            other._columns = 0;

            other._pitch = 0;
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

    /// @brief Gets values of specific row of matrix.
    /// @param row The index of row.
    /// @return The pointer to the values of row, aligned to a cache line boundary.
    auto
    data(const size_t row) -> double *
    {
        errors::assertMsgCritical(row < _rows, std::string("SimdMatrix.data: Index of row is out of range"));

        return _data + row * _pitch;
    }

    /// @brief Gets values of specific row of matrix.
    /// @param row The index of row.
    /// @return The constant pointer to the values of row, aligned to a cache line
    /// boundary.
    auto
    data(const size_t row) const -> const double *
    {
        errors::assertMsgCritical(row < _rows, std::string("SimdMatrix.data: Index of row is out of range"));

        return _data + row * _pitch;
    }

    /// @brief Gets number of rows in matrix.
    /// @return The number of rows.
    auto
    number_of_rows() const -> size_t
    {
        return _rows;
    }

    /// @brief Gets number of columns in matrix.
    /// @return The number of columns.
    auto
    number_of_columns() const -> size_t
    {
        return _columns;
    }

    /// @brief Gets padded number of columns in a row of matrix.
    /// @return The padded number of columns.
    auto
    pitch() const -> size_t
    {
        return _pitch;
    }

    /// @brief Gets number of values in matrix, padding included.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        return _rows * _pitch;
    }

    /// @brief Gets memory required to store the values of matrix.
    /// @return The memory in bytes.
    auto
    memory_size() const -> size_t
    {
        return number_of_elements() * sizeof(double);
    }

   private:
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

    /// @brief The values of matrix, stored row wise with padded rows.
    double *_data;

    /// @brief The number of rows in matrix.
    size_t _rows;

    /// @brief The number of columns in matrix.
    size_t _columns;

    /// @brief The padded number of columns in a row of matrix.
    size_t _pitch;
};

#endif /* SimdMatrix_hpp */
