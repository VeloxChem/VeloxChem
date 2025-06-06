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

#ifndef SubMatrix_hpp
#define SubMatrix_hpp

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

/// @brief Class CSubMatrix stores submatrix data and provides set of methods
/// for handling of submatrix data.
class CSubMatrix
{
   public:
    /// @brief The default constructor..
    CSubMatrix();

    /// @brief The constructor with submatrix dimensions.
    /// @param dimensions The dimensions of submatrix.
    CSubMatrix(const std::array<size_t, 4> &dimensions);

    /// @brief The constructor with submatrix dimensions and intial value.
    /// @param dimensions The dimensions of submatrix.
    /// @param factor The factor to set all values of submatrix.
    CSubMatrix(const std::array<size_t, 4> &dimensions, const double factor);

    /// @brief The constructor with submatrix dimensions and vector of initial values.
    /// @param values The vector of values to initialize submatrix values.
    /// @param dimensions The dimensions of submatrix.
    CSubMatrix(const std::vector<double> &values, const std::array<size_t, 4> &dimensions);

    /// @brief The default copy constructor.
    /// @param other The submatrix to be copied.
    CSubMatrix(const CSubMatrix &other);

    /// @brief The default move constructor.
    /// @param other The submatrix to be moved.
    CSubMatrix(CSubMatrix &&other) noexcept;

    /// @brief The default destructor.
    ~CSubMatrix() = default;

    /// @brief The default copy assignment operator.
    /// @param other The submatrix to be copy assigned.
    /// @return The assigned submatrix.
    auto operator=(const CSubMatrix &other) -> CSubMatrix &;

    /// @brief The default move assignment operator.
    /// @param other The submatrix to be move assigned.
    /// @return The assigned submatrix.
    auto operator=(CSubMatrix &&other) noexcept -> CSubMatrix &;

    /// @brief The equality operator.
    /// @param other The submatrix to be compared.
    /// @return True if submatrices are equal, False otherwise.
    auto operator==(const CSubMatrix &other) const -> bool;

    /// @brief The non-equality operator.
    /// @param other The submatrix to be compared.
    /// @return True if submatrices are not equal, False otherwise.
    auto operator!=(const CSubMatrix &other) const -> bool;

    /// @brief The addition operator.
    /// @param other The submatrix to be added.
    /// @return The sum of two submatrices.
    auto operator+(const CSubMatrix &other) const -> CSubMatrix;

    /// @brief Gets reference to specific submatirx element using global indexing
    /// scheme.
    /// @param index The {row, column} index in global indexing scheme.
    /// @return The reference to submatrix element.
    inline auto
    operator[](const std::pair<size_t, size_t> &index) -> double &
    {
        return _values[(index.first - _dimensions[0]) * _dimensions[3] + index.second - _dimensions[1]];
    }

    /// @brief Gets constant reference to specific submatirx element using global
    /// indexing scheme.
    /// @param index The {row, column} index in global indexing scheme.
    /// @return The constant reference to submatrix element.
    inline auto
    operator[](const std::pair<size_t, size_t> &index) const -> const double &
    {
        return _values[(index.first - _dimensions[0]) * _dimensions[3] + index.second - _dimensions[1]];
    }

    /// @brief Gets reference to specific submatirx element using local indexing
    /// scheme.
    /// @param index The {row, column} index in local indexing scheme.
    /// @return The reference to submatrix element.
    inline auto
    at(const std::pair<size_t, size_t> &index) -> double &
    {
        return _values[index.first * _dimensions[3] + index.second];
    }

    /// @brief Gets constant reference to specific submatirx element using local
    /// indexing scheme.
    /// @param index The {row, column} index in local indexing scheme.
    /// @return The constant reference to submatrix element.
    inline auto
    at(const std::pair<size_t, size_t> &index) const -> const double &
    {
        return _values[index.first * _dimensions[3] + index.second];
    }

    /// @brief Set offsets of submatrix in global indexing scheme.
    /// @param offsets The rows and collumns offsets in  {rows, columns} form.
    auto set_offsets(const std::pair<size_t, size_t> &offsets) -> void;

    /// Set values of submatrix.
    /// - Parameter values: the vector with submatrix values.

    /// @brief Set values of submatrix elements.
    /// @param values The vector of values to initialize submatrix elements.
    auto set_values(const std::vector<double> &values) -> void;

    /// @brief Set values of submatrix elements to zero.
    auto zero() -> void;

    /// @brief Scales values of submatrix elements by given factor.
    /// @param factor The scaling factor.
    auto scale(const double factor) -> void;

    /// Symmetrizes values of square submatrix (NOTE: stores 2 * value for
    /// diagonal).

    /// @brief Symmetrizes values of square submatrix.
    /// NOTE: diagonal elements 2 * value_ii and off-diagonal elements (value_ij
    /// + value_ji)
    auto symmetrize() -> void;

    /// @brief Getter for dimensions of submatrix.
    /// @return The dimensions of submatrix.
    auto get_dimensions() const -> std::array<size_t, 4>;

    /// @brief Gets values of all submatrix elements as flat vector.
    /// @return The vector with submatrix elements.
    auto get_values() const -> std::vector<double>;

    /// @brief Gets offset of matrix rows in global indexing scheme.
    /// @return The offset of matrix rows.
    inline auto
    offset_of_rows() const -> size_t
    {
        return _dimensions[0];
    }

    /// @brief Gets offset of matrix columns in global indexing scheme.
    /// @return  The offset of matrix columns.
    inline auto
    offset_of_columns() const -> size_t
    {
        return _dimensions[1];
    }

    /// @brief Gets number of rows in submatrix.
    /// @return The number of rows.
    inline auto
    number_of_rows() const -> size_t
    {
        return _dimensions[2];
    }

    /// @brief Gets number of columns in submatrix.
    /// @return The number of columns.
    inline auto
    number_of_columns() const -> size_t
    {
        return _dimensions[3];
    }

    /// @brief Gets number of elements in submatrix.
    /// @return The number of elements.
    inline auto
    number_of_elements() const -> size_t
    {
        return _dimensions[2] * _dimensions[3];
    }

    /// @brief Checks if submatrix is square matrix.
    /// @return True if submatrix is square matrix, False otherwise.
    inline auto
    is_square() const -> bool
    {
        return _dimensions[2] == _dimensions[3];
    }
    
    /// @brief Gets constant pointer to raw submatrix data.
    /// @return The constant pointer to raw submatrix data.
    inline auto
    data() const -> const double*
    {
        return _values.data();
    }
    
    /// @brief Gets pointer to raw submatrix data.
    /// @return The pointer to raw submatrix data.
    inline auto
    data()  -> double*
    {
        return _values.data();
    }

   private:
    /// @brief The dimensions of submatrix: row and column offsets and
    /// dimensions.
    std::array<size_t, 4> _dimensions;

    /// @brief The flat vector storage for submatrix values.
    std::vector<double> _values;
};

#endif /* SubMatrix_hpp */
