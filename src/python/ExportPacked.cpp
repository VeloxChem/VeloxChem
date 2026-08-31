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



#include "ExportPacked.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "ExportGeneral.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "PackedMatrix.hpp"

namespace vlx_packed {  // vlx_packed namespace

/// @brief Expands a packed matrix into the dense matrix as a numpy array.
/// @param matrix The packed matrix to expand.
/// @param max_memory The maximum memory of the dense matrix in gigabytes.
/// @return The dense matrix.
static auto
dense_to_numpy(const CPackedMatrix &matrix, const double max_memory) -> py::array_t<double>
{
    const auto nrows = matrix.number_of_rows();

    const auto ncols = matrix.number_of_columns();

    // NOTE: a triangular matrix is expanded into both of its triangles, so the
    // dense matrix is up to twice the matrix and its memory is checked against a
    // limit before it is allocated. The limit is enforced with an exception
    // rather than with a critical error, as the latter would terminate the
    // Python interpreter.

    const auto memory = static_cast<double>(nrows) * static_cast<double>(ncols) * static_cast<double>(sizeof(double));

    if (const auto limit = max_memory * 1024.0 * 1024.0 * 1024.0; memory > limit)
    {
        throw std::runtime_error("PackedMatrix.to_numpy: Dense matrix of " + std::to_string(memory / (1024.0 * 1024.0 * 1024.0)) +
                                 " GB exceeds limit of " + std::to_string(max_memory) + " GB");
    }

    auto array = py::array_t<double>(std::vector<py::ssize_t>{static_cast<py::ssize_t>(nrows), static_cast<py::ssize_t>(ncols)});

    matrix.to_dense(array.mutable_data());

    return array;
}

/// @brief Expands a packed matrix into the dense matrix in an array given by the
/// caller.
/// @param matrix The packed matrix to expand.
/// @param array The array to expand the matrix into.
/// @note The array is taken without conversion, as an array of another type
/// would be converted into a temporary and the expansion would fill that rather
/// than the array of the caller, leaving it untouched and saying nothing.
/// @note The dense matrix of a large basis is several gigabytes, which the
/// allocator returns to the system when it is freed, so a caller expanding one
/// matrix after another pays for faulting in every page of it every time.
/// Filling an array which the caller keeps avoids that, and is the reason this
/// exists beside to_numpy.
static auto
dense_to_numpy(const CPackedMatrix &matrix, py::array_t<double> &array) -> void
{
    const auto nrows = matrix.number_of_rows();

    const auto ncols = matrix.number_of_columns();

    // NOTE: the array is checked here rather than by a critical error, as the
    // latter would terminate the Python interpreter, and the expansion would
    // write outside the array of a caller which got the shape wrong.

    if ((array.ndim() != 2) || (static_cast<size_t>(array.shape(0)) != nrows) || (static_cast<size_t>(array.shape(1)) != ncols))
    {
        throw std::runtime_error("PackedMatrix.fill_numpy: Array must have " + std::to_string(nrows) + " rows and " +
                                 std::to_string(ncols) + " columns");
    }

    if ((array.flags() & py::array::c_style) == 0)
    {
        throw std::runtime_error("PackedMatrix.fill_numpy: Array must be C contiguous");
    }

    matrix.to_dense(array.mutable_data());
}

auto
export_packed(py::module &m) -> void
{
    // CPackedMatrix class

    PyClass<CPackedMatrix>(m, "PackedMatrix")
        .def(py::init<>())
        .def(py::init<const size_t, const size_t, const mat_t>(),
             "Creates a packed matrix of given dimensions.",
             py::arg("nrows"),
             py::arg("ncols"),
             py::arg("matrix_type"))
        .def(py::init<const CMolecularBasis &, const mat_t>(),
             "Creates a packed matrix for a molecular basis.",
             py::arg("basis"),
             py::arg("matrix_type"))
        .def(py::init<const CMolecularBasis &, const CMolecularBasis &>(),
             "Creates a general packed matrix for a pair of molecular bases.",
             py::arg("bra_basis"),
             py::arg("ket_basis"))
        .def("number_of_rows", &CPackedMatrix::number_of_rows, "Gets number of rows of matrix.")
        .def("number_of_columns", &CPackedMatrix::number_of_columns, "Gets number of columns of matrix.")
        .def("get_type", &CPackedMatrix::get_type, "Gets type of matrix.")
        .def("is_triangular", &CPackedMatrix::is_triangular, "Checks if only one triangle of matrix is stored.")
        .def("number_of_elements", &CPackedMatrix::number_of_elements, "Gets number of values stored by matrix.")
        .def("memory_size", &CPackedMatrix::memory_size, "Gets memory required to store the values of matrix in bytes.")
        .def("zero", &CPackedMatrix::zero, "Sets the values of matrix to zero.")
        .def("scale", &CPackedMatrix::scale, "Scales the values of matrix.", py::arg("factor"))
        .def("index",
             &CPackedMatrix::index,
             "Gets position of the value storing the element of matrix.",
             py::arg("irow"),
             py::arg("icol"))
        .def("at", &CPackedMatrix::at, "Gets element of matrix.", py::arg("irow"), py::arg("icol"))
        .def(
            "values_to_numpy",
            [](const CPackedMatrix &self) -> py::array_t<double> {
                return py::array_t<double>(std::vector<py::ssize_t>{static_cast<py::ssize_t>(self.number_of_elements())}, self.data());
            },
            "Gets the values of matrix as they are stored, i.e. the lower triangle of a triangular matrix.")
        .def(
            "to_numpy",
            [](const CPackedMatrix &self, const double max_memory) -> py::array_t<double> {
                return vlx_packed::dense_to_numpy(self, max_memory);
            },
            "Gets the dense matrix in atomic orbital basis. The matrix may be large.",
            py::arg("max_memory") = 8.0)
        .def(
            "fill_numpy",
            [](const CPackedMatrix &self, py::array_t<double> &array) -> void { vlx_packed::dense_to_numpy(self, array); },
            "Expands the matrix into the dense matrix in atomic orbital basis in the array of the caller.",
            py::arg("array").noconvert());
}

}  // namespace vlx_packed
