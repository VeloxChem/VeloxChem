//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "ExportMath.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "AngularMomentum.hpp"
#include "CantorFunc.hpp"
#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MatrixIndex.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_math {  // vlx_math namespace

// Helper function for CDenseMatrix constructor

auto
CDenseMatrix_from_numpy(const py::array_t<double>& arr) -> std::shared_ptr<CDenseMatrix>
{
    // check dimension

    std::string errdim("DenseMatrix: Expecting a 2D numpy array");

    errors::assertMsgCritical(arr.ndim() == 2, errdim);

    if (arr.data() == nullptr || arr.size() == 0)
    {
        return std::make_shared<CDenseMatrix>();
    }

    // check that the numpy array is c-style contiguous

    std::string errsrc("DenseMatrix: Expecting a contiguous numpy array");

    auto c_style = py::detail::check_flags(arr.ptr(), py::array::c_style);

    auto f_style = py::detail::check_flags(arr.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style | f_style, errsrc);

    // create CDenseMatrix from numpy array

    int64_t nrows = static_cast<int64_t>(arr.shape(0));

    int64_t ncols = static_cast<int64_t>(arr.shape(1));

    T4Index dim({0, 0, nrows, ncols});

    CSubMatrix submat(dim);

    auto submat_ptr = submat.getData();

    if (c_style)
    {
        std::memcpy(submat_ptr, arr.data(), arr.size() * sizeof(double));
    }
    else if (f_style)
    {
        for (py::ssize_t i = 0; i < arr.shape(0); i++)
        {
            for (py::ssize_t j = 0; j < arr.shape(1); j++)
            {
                submat_ptr[i * arr.shape(1) + j] = arr.data()[j * arr.shape(0) + i];
            }
        }
    }

    return std::make_shared<CDenseMatrix>(submat);
}

// Exports classes/functions in src/math to python

auto
export_math(py::module& m) -> void
{
    // exposing enum from MatrixType.hpp

    // clang-format off
    py::enum_<mat_t>(m, "mat_t")
        .value("symm", mat_t::symm)
        .value("antisymm", mat_t::antisymm)
        .value("gen", mat_t::gen);
    // clang-format on

    // exposing functions from MathConst.hpp

    m.def("get_pi", &mathconst::getPiValue, "Gets PI value.");

    // exposing functions from MathIndex.hpp

    m.def("uplo_index", &mathfunc::uplo_index, "Gets index of upper triangular matrix.");

    // exposing functions from CantorFunc.hpp

    m.def("cantor_index", &mathfunc::getCantorIndex, "Computes Cantor index for pair of non-negative numbers.");

    m.def("cantor_pair", &mathfunc::getCantorPair, "Extractes pair of non-negative numbers from Cantor index.");

    // exposing functions from AngularMomentum.hpp

    m.def("to_spherical_components",
          py::overload_cast<int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given angular momentum.");

    m.def("to_spherical_components",
          py::overload_cast<int64_t, int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given pair of angular momentums.");

    m.def("to_cartesian_components",
          py::overload_cast<int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given angular momentum.");

    m.def("to_cartesian_components",
          py::overload_cast<int64_t, int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given pair of angular momentums.");

    m.def("angular_component_to_str", &angmom::getStringOfAngularMomentum, "Gets string of angular momentum component.");

    // exposing functions from MatrixFunc.hpp

    m.def(
        "make_matrix",
        [](const CMolecularBasis& basis, const mat_t mat_type) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::makeMatrix(basis, mat_type));
        },
        "Creates matrix for given basis.");

    m.def(
        "make_matrix",
        [](const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::makeMatrix(bra_basis, ket_basis));
        },
        "Creates matrix for given pair of bases.");

    // CSubMatrix class

    PyClass<CSubMatrix>(m, "SubMatrix")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def(py::init<const T4Index&>())
        .def(py::init<const std::vector<double>&, const T4Index&>())
        .def(
            "at",
            [](CSubMatrix& self, const int64_t irow, const int64_t icol, const double value) -> void { self.at(irow, icol, true) = value; },
            "Sets value of submatrix element using supermatrix indexing scheme.")
        .def(
            "at",
            [](const CSubMatrix& self, const int64_t irow, const int64_t icol) -> double { return self.at(irow, icol, true); },
            "Gets value of submatrix element using supermatrix indexing scheme.")
        .def(
            "set_value",
            [](CSubMatrix& self, const int64_t irow, const int64_t icol, const double value) -> void { self.at(irow, icol, false) = value; },
            "Sets value of submatrix element.")
        .def(
            "get_value",
            [](const CSubMatrix& self, const int64_t irow, const int64_t icol) -> double { return self.at(irow, icol, false); },
            "Gets value of submatrix element.")
        .def("set_offsets", &CSubMatrix::setOffsets, "Sets offsets for rows and columns of submatrix in supermatrix.")
        .def(
            "set_values",
            [](CSubMatrix& self, const py::array_t<double>& values) -> void {
                if (values.ndim() == 2)
                {
                    const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());

                    const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());

                    if ((nrows == values.shape(0)) && (ncols == values.shape(1)))
                    {
                        std::memcpy(self.getData(), values.data(), nrows * ncols * sizeof(double));
                    }
                }
            },
            "Sets values of submatrix using numpy array.")
        .def("zero", &CSubMatrix::zero, "Sets values of submatrix to zero.")
        .def("symmetrize", &CSubMatrix::symmetrize, "Symmetrizes values of square submatrix.")
        .def(
            "to_numpy",
            [](const CSubMatrix& self) -> py::array_t<double> {
                const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());

                const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());

                const auto tdim = static_cast<py::ssize_t>(sizeof(double));

                return py::array_t<double>(std::vector<py::ssize_t>({nrows, ncols}), std::vector<py::ssize_t>({ncols * tdim, tdim}), self.getData());
            },
            "Converts submatrix to numpy array.")
        .def("get_dimensions", &CSubMatrix::getDimensions, "Gets dimensions of submatrix.")
        .def("offset_of_rows", &CSubMatrix::getOffsetOfRows, "Gets offset of rows in submatrix.")
        .def("offset_of_columns", &CSubMatrix::getOffsetOfColumns, "Gets offset of columns in submatrix.")
        .def("number_of_rows", &CSubMatrix::getNumberOfRows, "Gets number of rows in submatrix.")
        .def("number_of_columns", &CSubMatrix::getNumberOfColumns, "Gets number of columns in submatrix.");

    // CMatrix class

    PyClass<CMatrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<const std::map<T2Pair, CSubMatrix>&, const mat_t>())
        .def(py::init<const CMatrix&>())
        .def("add", py::overload_cast<const CSubMatrix&, const T2Pair&>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("add", py::overload_cast<const T4Index&, const T2Pair&>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("set_type", &CMatrix::setType, "Sets matrix type.")
        .def("zero", &CMatrix::zero, "Sets values of matrix to zero.")
        .def(
            "set_values",
            [](CMatrix& self, const py::array_t<double>& values) -> void {
                if (values.ndim() == 2)
                {
                    const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());

                    const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());

                    if ((nrows == values.shape(0)) && (ncols == values.shape(1)))
                    {
                        for (auto tpair : self.getAngularPairs())
                        {
                            auto submat = self.getSubMatrix(tpair);

                            const auto [row_off, col_off, nrows, ncols] = submat->getDimensions();

                            for (int64_t i = 0; i < nrows; i++)
                            {
                                for (int64_t j = 0; j < ncols; j++)
                                {
                                    submat->at(i, j, false) = values.at(i + row_off, j + col_off);
                                }
                            }
                        }
                    }
                }
            },
            "Sets values of matrix using numpy array.")
        .def("get_type", &CMatrix::getType, "Gets matrix type.")
        .def("get_angular_pairs", &CMatrix::getAngularPairs, "Gets vector of angular pairs for stored submatrices.")
        .def(
            "get_submatrix",
            [](const CMatrix& self, const T2Pair& angpair) -> std::shared_ptr<CSubMatrix> {
                if (auto submat = self.getSubMatrix(angpair); submat != nullptr)
                {
                    return std::make_shared<CSubMatrix>(*submat);
                }
                else
                {
                    return std::make_shared<CSubMatrix>();
                }
            },
            "Gets specific submatrix from matrix.")
        .def("is_angular_order", &CMatrix::isAngularOrder, "Checks if submatrix with this angular pair is stored in matrix.")
        .def("number_of_rows", &CMatrix::getNumberOfRows, "Number of rows in matrix.")
        .def("number_of_columns", &CMatrix::getNumberOfColumns, "Number of columns in matrix.")
        .def(
            "get_full_matrix",
            [](const CMatrix& self) -> std::shared_ptr<CSubMatrix> { return std::make_shared<CSubMatrix>(self.getFullMatrix()); },
            "Creates full matrix representation of matrix.");

    // CDenseMatrix class

    PyClass<CDenseMatrix>(m, "DenseMatrix")
        .def(py::init<>())
        .def(py::init<const int64_t, const int64_t>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init(&CDenseMatrix_from_numpy))
        .def("number_of_rows", &CDenseMatrix::getNumberOfRows, "Gets number of rows in dense matrix.")
        .def("number_of_columns", &CDenseMatrix::getNumberOfColumns, "Gets number of columns in dense matrix.")
        .def("symmetrize", &CDenseMatrix::symmetrize, "Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).")
        .def("slice",
             &CDenseMatrix::slice,
             "Creates dense matrix object by slicing columns at selected position from this dense matrix object.",
             "i_column"_a,
             "n_columns"_a)
        .def(
            "to_numpy",
            [](const CDenseMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.values(), std::vector<int64_t>{self.getNumberOfRows(), self.getNumberOfColumns()});
            },
            "Converts DenseMatrix to numpy array.")
        .def(py::self == py::self);
}

}  // namespace vlx_math
