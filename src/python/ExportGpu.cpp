//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "ExportGpu.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "CudaDevices.hpp"
#include "ErrorHandler.hpp"
#include "FockDriverGPU.hpp"
#include "ExportGeneral.hpp"
#include "ScreeningData.hpp"
#include "XCIntegratorGPU.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_gpu {  // vlx_gpu namespace

// Exports classes/functions in src/gpu to python

void
export_gpu(py::module& m)
{
    // CCudaDevices class

    py::class_<CCudaDevices, std::shared_ptr<CCudaDevices>>(m, "CudaDevices")
        .def(py::init<>())
        .def("get_number_devices", &CCudaDevices::getNumberOfDevices)
        .def("__str__", &CCudaDevices::getString);

    // CScreeningData class

    py::class_<CScreeningData, std::shared_ptr<CScreeningData>>(m, "ScreeningData")
        .def(py::init<const CMolecule&, const CMolecularBasis&>());

    m.def(
        "compute_gto_values",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> py::array_t<double> {
            auto gtovalues = gpu::computeGtoValuesOnGridPoints(molecule, basis, molecularGrid);
            return vlx_general::pointer_to_numpy(gtovalues.values(), {gtovalues.getNumberOfRows(), gtovalues.getNumberOfColumns()});
        },
        "Computes GTO values on grid points using GPU.");

    m.def(
        "compute_gto_values_and_derivatives",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> py::list {
            auto     gto_values_derivs = gpu::computeGtoValuesAndDerivativesOnGridPoints(molecule, basis, molecularGrid);
            py::list ret;
            for (size_t i = 0; i < gto_values_derivs.size(); i++)
            {
                ret.append(vlx_general::pointer_to_numpy(gto_values_derivs[i].values(),
                                                         {gto_values_derivs[i].getNumberOfRows(), gto_values_derivs[i].getNumberOfColumns()}));
            }
            return ret;
        },
        "Computes GTO values and derivatives on grid points using GPU.");

    m.def(
        "add_matrix_gpu",
        [](py::array_t<double>& Y, const std::vector<py::array_t<double>>& Xs, const std::vector<double>& alphas) -> void {
            //std::string errshape("add_matrix_gpu: Mismatch in shape");
            //std::string errstyle("add_matrix_gpu: Expecting contiguous numpy array");

            //auto c_style_Y = py::detail::check_flags(Y.ptr(), py::array::c_style);
            //errors::assertMsgCritical(c_style_Y, errstyle);

            std::vector<const double*> X_ptr_s;
            for (auto& X : Xs)
            {
                //auto c_style_X = py::detail::check_flags(X.ptr(), py::array::c_style);
                //errors::assertMsgCritical(c_style_X, errstyle);
                //errors::assertMsgCritical(Y.size() == X.size(), errshape);
                X_ptr_s.push_back(X.data());
            }

            const auto n2 = static_cast<int64_t>(Y.size());
            gpu::addMatrix(Y.mutable_data(), X_ptr_s, n2, alphas);
        },
        "Adds matrix using GPU.");

    m.def(
        "matmul_gpu",
        [](const py::array_t<double>& A, const py::array_t<double>& B) -> py::array_t<double> {
            std::string errshape("matmul_gpu: Mismatch in matrix shape");

            errors::assertMsgCritical(A.shape(1) == B.shape(0), errshape);

            std::string errstyle("matmul_gpu: Expecting contiguous numpy array");

            auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);
            auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);

            errors::assertMsgCritical(c_style_A && c_style_B, errstyle);

            const auto nrows_A = static_cast<int64_t>(A.shape(0));
            const auto ncols_A = static_cast<int64_t>(A.shape(1));
            const auto ncols_B = static_cast<int64_t>(B.shape(1));

            std::vector<double> C(nrows_A * ncols_B);

            gpu::computeMatrixMultiplication(C.data(), A.data(), B.data(), nrows_A, ncols_A, ncols_B);

            return vlx_general::pointer_to_numpy(C.data(), {nrows_A, ncols_B});
        },
        "Computes matrix multiplication using GPU.");

    m.def(
        "eigh_gpu",
        [](const py::array_t<double>& A) -> py::list {
            // check dimension and shape

            errors::assertMsgCritical(A.ndim() == 2, "eigh_gpu: Invalid shape of matrix A");

            auto nrows_A = A.shape(0);
            auto ncols_A = A.shape(1);

            errors::assertMsgCritical(ncols_A == nrows_A, "eigh_gpu: Matrix A is not symmetric");

            // check layout

            auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

            errors::assertMsgCritical(c_style_A, "eigh_gpu: Matrix A is not C-style contiguous");

            // initialize eigenvalues and eigenvectors

            auto dim = nrows_A;

            py::array_t<double> eigenValues(dim);
            py::array_t<double> eigenVectors({dim, dim});

            auto evecs = eigenVectors.mutable_data();
            auto evals = eigenValues.mutable_data();

            std::memcpy(evecs, A.data(), A.size() * sizeof(double));

            // diagonalize matrix

            gpu::diagonalizeMatrix(evecs, evals, static_cast<int64_t>(nrows_A));

            py::list result;

            result.append(eigenValues);

            result.append(eigenVectors);

            return result;

            },
        "Diagonalizes matrix using GPU.");

    m.def("integrate_vxc_fock", &gpu::integrateVxcFock, "Integrates Vxc matrix using GPU.");

    m.def("compute_fock_gpu", &gpu::computeFockOnGPU, "Computes Fock matrix using GPU.");

    m.def("compute_one_electron_integrals_gpu", &gpu::computeOneElectronIntegralsOnGPU, "Computes one-electron integral matrices using GPU.");
}

}  // namespace vlx_gpu
