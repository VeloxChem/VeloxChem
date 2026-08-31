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



#ifndef PackedMatrix_hpp
#define PackedMatrix_hpp

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"

/// @brief Class CPackedMatrix stores the values of a matrix in atomic orbital
/// basis in full, without a sparsity pattern. It is the storage of the integrals
/// which cannot be screened, such as the two-center electron repulsion
/// integrals, whose slow decay with the interatomic distance leaves no atom pair
/// below the threshold and whose matrix is dense for every molecule.
/// @note A symmetric or antisymmetric matrix is stored as its lower triangle
/// including the diagonal, in row major order, so the element of row i and
/// column j with j <= i is at i (i + 1) / 2 + j. The offset of a row does not
/// depend on the dimensions of the matrix and each row is contiguous. The
/// diagonal of an antisymmetric matrix is stored although it is zero, so that
/// the two types share one layout and one index and differ only in the sign
/// carried by the elements above the diagonal.
/// @note A general matrix is stored in full, in row major order, as neither of
/// its triangles determines the other.
/// @note The values are held in one allocation, which the constructor makes and
/// zeroes, unlike CSparseMatrix, whose values are allocated one block at a time
/// and after its sparsity patterns are set up. There is one allocation here and
/// nothing to describe before it is made, so deferring it would buy nothing.
class CPackedMatrix
{
   public:
    /// @brief The number of values one thread zeroes, scales or copies at a
    /// time. The matrix reaches several gigabytes, so these traversals are
    /// divided over the threads rather than left to the calling one.
    static constexpr size_t _chunk_size = 1 << 20;

    /// @brief The default constructor.
    CPackedMatrix()

        : _nrows(0)

        , _ncols(0)

        , _type(mat_t::general)

        , _values{}
    {
    }

    /// @brief The constructor with dimensions and matrix type.
    /// @param nrows The number of rows of matrix.
    /// @param ncols The number of columns of matrix.
    /// @param mat_type The type of matrix.
    CPackedMatrix(const size_t nrows, const size_t ncols, const mat_t mat_type)

        : _nrows(nrows)

        , _ncols(ncols)

        , _type(mat_type)

        , _values{}
    {
        if ((mat_type != mat_t::general) && (nrows != ncols))
        {
            errors::assertMsgCritical(false, std::string("PackedMatrix: Symmetric and antisymmetric matrix must be square"));
        }

        _values.resize(_number_of_elements(nrows, ncols, mat_type), 0.0);
    }

    /// @brief The constructor with molecular basis and matrix type.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param mat_type The type of matrix.
    CPackedMatrix(const CMolecularBasis &basis, const mat_t mat_type)

        : CPackedMatrix(basis.dimensions_of_basis(), basis.dimensions_of_basis(), mat_type)
    {
    }

    /// @brief The constructor with a pair of molecular bases. The matrix of two
    /// molecular bases is general, as its rows and columns are indexed by
    /// different basis functions and neither of its triangles determines the
    /// other.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    CPackedMatrix(const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis)

        : CPackedMatrix(bra_basis.dimensions_of_basis(), ket_basis.dimensions_of_basis(), mat_t::general)
    {
    }

    /// @brief Gets number of rows of matrix.
    /// @return The number of rows.
    auto
    number_of_rows() const -> size_t
    {
        return _nrows;
    }

    /// @brief Gets number of columns of matrix.
    /// @return The number of columns.
    auto
    number_of_columns() const -> size_t
    {
        return _ncols;
    }

    /// @brief Gets type of matrix.
    /// @return The type of matrix.
    auto
    get_type() const -> mat_t
    {
        return _type;
    }

    /// @brief Checks if only one triangle of matrix is stored.
    /// @return True if the matrix is stored as its lower triangle, false if it
    /// is stored in full.
    auto
    is_triangular() const -> bool
    {
        return _type != mat_t::general;
    }

    /// @brief Gets number of values stored by matrix.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        return _values.size();
    }

    /// @brief Gets memory required to store the values of matrix.
    /// @return The memory in bytes.
    auto
    memory_size() const -> size_t
    {
        return _values.size() * sizeof(double);
    }

    /// @brief Gets values of matrix.
    /// @return The pointer to the values.
    auto
    data() -> double *
    {
        return _values.data();
    }

    /// @brief Gets values of matrix.
    /// @return The constant pointer to the values.
    auto
    data() const -> const double *
    {
        return _values.data();
    }

    /// @brief Gets position of the value storing the element of matrix.
    /// @param irow The index of row of the element.
    /// @param icol The index of column of the element.
    /// @return The position of the value among the values of matrix.
    /// @note The element above the diagonal of a triangular matrix is stored by
    /// the value of the transposed element, so the position is the same for both
    /// of them and the sign which distinguishes them is not applied here.
    /// @note The indices are not checked, as this is called once per element of
    /// a matrix which has no sparsity pattern to make it smaller. Use at() to
    /// read a single element.
    auto
    index(const size_t irow, const size_t icol) const -> size_t
    {
        if (_type == mat_t::general) return irow * _ncols + icol;

        const auto i = (irow < icol) ? icol : irow;

        const auto j = (irow < icol) ? irow : icol;

        return i * (i + 1) / 2 + j;
    }

    /// @brief Gets element of matrix.
    /// @param irow The index of row of the element.
    /// @param icol The index of column of the element.
    /// @return The value of the element.
    /// @note The element above the diagonal of a triangular matrix is not
    /// stored, so it is taken from the transposed element with the sign of the
    /// matrix type.
    auto
    at(const size_t irow, const size_t icol) const -> double
    {
        if ((irow >= _nrows) || (icol >= _ncols))
        {
            errors::assertMsgCritical(false, std::string("PackedMatrix.at: Index of element is out of range"));
        }

        const auto fval = _values[index(irow, icol)];

        return ((_type == mat_t::antisymmetric) && (irow < icol)) ? -fval : fval;
    }

    /// @brief Sets the values of matrix to zero.
    auto
    zero() -> void
    {
        const auto nchunks = _number_of_chunks();

#pragma omp parallel for schedule(static) if (nchunks > 1)
        for (int i = 0; i < nchunks; i++)
        {
            const auto [first, last] = _chunk_range(i);

            std::fill(_values.data() + first, _values.data() + last, 0.0);
        }
    }

    /// @brief Scales the values of matrix.
    /// @param factor The scaling factor.
    auto
    scale(const double factor) -> void
    {
        const auto nchunks = _number_of_chunks();

#pragma omp parallel for schedule(static) if (nchunks > 1)
        for (int i = 0; i < nchunks; i++)
        {
            const auto [first, last] = _chunk_range(i);

            for (size_t k = first; k < last; k++)
            {
                _values[k] *= factor;
            }
        }
    }

    /// @brief Expands matrix into the dense matrix in atomic orbital basis.
    /// @param values The values of the dense matrix, as a row major array of
    /// number_of_rows() rows and number_of_columns() columns, which must be
    /// allocated by the caller.
    /// @note The elements above the diagonal of a triangular matrix are the
    /// transposed elements, with their sign for an antisymmetric matrix.
    auto
    to_dense(double *values) const -> void
    {
        if (_type == mat_t::general)
        {
            const auto nchunks = _number_of_chunks();

#pragma omp parallel for schedule(static) if (nchunks > 1)
            for (int i = 0; i < nchunks; i++)
            {
                const auto [first, last] = _chunk_range(i);

                std::copy(_values.data() + first, _values.data() + last, values + first);
            }

            return;
        }

        const auto fsign = (_type == mat_t::antisymmetric) ? -1.0 : 1.0;

        const auto nrows = static_cast<int>(_nrows);

        // NOTE: the rows are interleaved over the threads, as the elements of a
        // row up to the diagonal are copied from the contiguous values of that
        // row while those beyond it are gathered from the column of the same
        // index. The first rows are therefore gathered almost entirely and the
        // last are copied almost entirely, and a thread given a contiguous
        // stretch of the rows would get either the one or the other.

#pragma omp parallel for schedule(static, 1) if (nrows > 1)
        for (int i = 0; i < nrows; i++)
        {
            const auto irow = static_cast<size_t>(i);

            auto *row = values + irow * _ncols;

            const auto *packed = _values.data() + irow * (irow + 1) / 2;

            std::copy(packed, packed + irow + 1, row);

            for (size_t j = irow + 1; j < _ncols; j++)
            {
                row[j] = fsign * _values[j * (j + 1) / 2 + irow];
            }
        }
    }

   private:
    /// @brief Gets number of values required to store a matrix.
    /// @param nrows The number of rows of matrix.
    /// @param ncols The number of columns of matrix.
    /// @param mat_type The type of matrix.
    /// @return The number of values.
    static auto
    _number_of_elements(const size_t nrows, const size_t ncols, const mat_t mat_type) -> size_t
    {
        return (mat_type == mat_t::general) ? nrows * ncols : nrows * (nrows + 1) / 2;
    }

    /// @brief Gets number of chunks the values of matrix are traversed in.
    /// @return The number of chunks.
    auto
    _number_of_chunks() const -> int
    {
        return static_cast<int>((_values.size() + _chunk_size - 1) / _chunk_size);
    }

    /// @brief Gets range of the values of matrix covered by a chunk.
    /// @param index The index of chunk.
    /// @return The first value of the chunk and the value past its last.
    auto
    _chunk_range(const int index) const -> std::pair<size_t, size_t>
    {
        const auto first = static_cast<size_t>(index) * _chunk_size;

        return {first, std::min(first + _chunk_size, _values.size())};
    }

    /// @brief The number of rows of matrix.
    size_t _nrows;

    /// @brief The number of columns of matrix.
    size_t _ncols;

    /// @brief The type of matrix.
    mat_t _type;

    /// @brief The values of matrix.
    std::vector<double> _values;
};

#endif /* PackedMatrix_hpp */
