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

#include "ExportGpu.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "FockDriverGPU.hpp"
#include "FockGradientDriverGPU.hpp"
#include "GpuDevices.hpp"
#include "LinearAlgebraGPU.hpp"
#include "ScreeningData.hpp"
#include "GradientScreeningData.hpp"
#include "XCIntegratorGPU.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_gpu {  // vlx_gpu namespace

// Exports classes/functions in src/gpu to python

void
export_gpu(py::module& m)
{
    // CGpuDevices class

    py::class_<CGpuDevices, std::shared_ptr<CGpuDevices>>(m, "GpuDevices")
        .def(py::init<>())
        .def("get_number_devices", &CGpuDevices::getNumberOfDevices)
        .def("__str__", &CGpuDevices::getString);

    // CScreeningData class

    py::class_<CScreeningData, std::shared_ptr<CScreeningData>>(m, "ScreeningData")
        .def(py::init<const CMolecule&, const CMolecularBasis&, const int64_t, const double, const double, const int, const int>())
        .def("get_num_gpus_per_node", &CScreeningData::getNumGpusPerNode)
        .def("reset_mpi", &CScreeningData::reset_mpi)
        .def("get_timer_summary", &CScreeningData::getTimerSummary)
        .def("get_gpu_timer_summary", &CScreeningData::getGpuTimerSummary)
        .def("get_prim_cart_naos", &CScreeningData::get_number_of_prim_cart_aos)
        .def(
            "fill_prim_q_and_cart_dens",
            [](CScreeningData&         self,
               const CMolecule&        mol,
               const CMolecularBasis&  bas,
               const CAODensityMatrix& dmat,
               CDenseMatrix&           prim_q_mat,
               CDenseMatrix&           prim_cart_d_mat) -> void {
                const auto cart_dmat = gpu::transformDensity(mol, bas, dmat);
                const auto cart_naos = cart_dmat.getNumberOfRows();
                const auto cart_dens_ptr = cart_dmat.values();
                return self.fill_mat_Q_and_D(cart_naos, cart_dens_ptr, prim_q_mat, prim_cart_d_mat);
            },
            "Gets primitive Cartesian density matrix.")
        .def("compute_q_prime_slice", &CScreeningData::get_Q_prime_slice);

    // CGradientScreeningData class

    py::class_<CGradientScreeningData, std::shared_ptr<CGradientScreeningData>>(m, "GradientScreeningData")
        .def(py::init<const CMolecule&, const CMolecularBasis&, const CAODensityMatrix&, const CDenseMatrix&, const int64_t, const double, const double, const int, const int>())
        .def("get_num_gpus_per_node", &CGradientScreeningData::getNumGpusPerNode);

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

    m.def("get_avail_mem_gpu", &gpu::getAvailableGpuMemory);

    m.def(
        "dot_product_gpu",
        [](const py::array_t<double>& A, const py::array_t<double>& B) -> double {
            std::string errsize("dot_product_gpu: Mismatch in size");
            errors::assertMsgCritical(A.size() == B.size(), errsize);

            std::string errstyle("dot_product_gpu: Expecting contiguous numpy array");
            auto        c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);
            auto        c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);
            errors::assertMsgCritical(c_style_A && c_style_B, errstyle);

            const auto n = static_cast<int64_t>(A.size());

            return gpu::computeDotProduct(A.data(), B.data(), n);
        },
        "Computes dot product.");

    m.def(
        "weighted_sum_gpu",
        [](const std::vector<double>& weights, const std::vector<py::array_t<double>>& arrays) -> py::array_t<double> {
            std::string errsize("weighted_sum_gpu: Mismatch in size");
            errors::assertMsgCritical(weights.size() == arrays.size(), errsize);

            std::string errzero("weighted_sum_gpu: Expecting at least one array");
            errors::assertMsgCritical(arrays.size() > 0, errzero);

            const auto nrows = static_cast<int64_t>(arrays[0].shape(0));
            const auto ncols = static_cast<int64_t>(arrays[0].shape(1));

            std::vector<double> weighted_array(nrows * ncols, 0.0);

            std::vector<const double*> data_pointers;

            for (size_t i = 0; i < arrays.size(); i++)
            {
                std::string errstyle("weighted_sum_gpu: Expecting contiguous numpy array");
                auto        c_style = py::detail::check_flags(arrays[i].ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);

                data_pointers.push_back(arrays[i].data());
            }

            gpu::computeWeightedSum(weighted_array.data(), weights, data_pointers, nrows * ncols);

            return vlx_general::pointer_to_numpy(weighted_array.data(), {nrows, ncols});
        },
        "Computes weighted sum of matrices.");

    m.def(
        "compute_error_vector_gpu",
        [](const py::array_t<double>& X, const py::array_t<double>& F, const py::array_t<double>& D, const py::array_t<double>& S)
            -> py::array_t<double> {
            std::string errshape("compute_error_vector_gpu: Mismatch in matrix shape");
            std::string errstyle("compute_error_vector_gpu: Expecting contiguous numpy array");

            errors::assertMsgCritical((X.shape(0) == F.shape(0)) && (X.shape(0) == F.shape(1)), errshape);
            errors::assertMsgCritical((X.shape(0) == D.shape(0)) && (X.shape(0) == D.shape(1)), errshape);
            errors::assertMsgCritical((X.shape(0) == S.shape(0)) && (X.shape(0) == S.shape(1)), errshape);

            auto c_style_X = py::detail::check_flags(X.ptr(), py::array::c_style);
            auto f_style_X = py::detail::check_flags(X.ptr(), py::array::f_style);
            errors::assertMsgCritical(c_style_X | f_style_X, errstyle);

            const auto trans_X = c_style_X ? std::string("N") : std::string("T");

            auto c_style_F = py::detail::check_flags(F.ptr(), py::array::c_style);
            auto c_style_D = py::detail::check_flags(D.ptr(), py::array::c_style);
            auto c_style_S = py::detail::check_flags(S.ptr(), py::array::c_style);
            errors::assertMsgCritical(c_style_F & c_style_D & c_style_S, errstyle);

            const auto nao = static_cast<int64_t>(X.shape(0));
            const auto nmo = static_cast<int64_t>(X.shape(1));

            std::vector<double> errvec(nmo * nmo);

            gpu::computeErrorVector(errvec.data(), X.data(), F.data(), D.data(), S.data(), nmo, nao, trans_X);

            return vlx_general::pointer_to_numpy(errvec.data(), {nmo, nmo});
        },
        "Computes error vector using GPU.");

    m.def(
        "transform_matrix_gpu",
        [](const py::array_t<double>& X, const py::array_t<double>& F) -> py::array_t<double> {
            std::string errshape("transform_matrix_gpu: Mismatch in matrix shape");
            std::string errstyle("transform_matrix_gpu: Expecting contiguous numpy array");

            errors::assertMsgCritical((X.shape(0) == F.shape(0)) && (X.shape(0) == F.shape(1)), errshape);

            auto c_style_X = py::detail::check_flags(X.ptr(), py::array::c_style);
            auto f_style_X = py::detail::check_flags(X.ptr(), py::array::f_style);
            errors::assertMsgCritical(c_style_X | f_style_X, errstyle);

            const auto trans_X = c_style_X ? std::string("N") : std::string("T");

            auto c_style_F = py::detail::check_flags(F.ptr(), py::array::c_style);
            errors::assertMsgCritical(c_style_F, errstyle);

            const auto nao = static_cast<int64_t>(X.shape(0));
            const auto nmo = static_cast<int64_t>(X.shape(1));

            std::vector<double> transformed_F(nmo * nmo);

            gpu::transformMatrix(transformed_F.data(), X.data(), F.data(), nmo, nao, trans_X);

            return vlx_general::pointer_to_numpy(transformed_F.data(), {nmo, nmo});
        },
        "Transform matrix using GPU.");

    m.def(
        "matmul_gpu",
        [](const py::array_t<double>& A, const py::array_t<double>& B) -> py::array_t<double> {
            std::string errshape("matmul_gpu: Mismatch in matrix shape");
            std::string errstyle("matmul_gpu: Expecting contiguous numpy array");

            errors::assertMsgCritical(A.shape(1) == B.shape(0), errshape);

            auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);
            auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

            auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);
            auto f_style_B = py::detail::check_flags(B.ptr(), py::array::f_style);

            errors::assertMsgCritical(c_style_A | f_style_A, errstyle);
            errors::assertMsgCritical(c_style_B | f_style_B, errstyle);

            const auto trans_A = c_style_A ? std::string("N") : std::string("T");
            const auto trans_B = c_style_B ? std::string("N") : std::string("T");

            const auto nrows_A = static_cast<int64_t>(A.shape(0));
            const auto ncols_A = static_cast<int64_t>(A.shape(1));
            const auto ncols_B = static_cast<int64_t>(B.shape(1));

            std::vector<double> C(nrows_A * ncols_B);

            gpu::computeMatrixMultiplication(C.data(), A.data(), B.data(), trans_A, trans_B, nrows_A, ncols_A, ncols_B);

            return vlx_general::pointer_to_numpy(C.data(), {nrows_A, ncols_B});
        },
        "Computes matrix multiplication using GPU.");

    m.def(
        "eigh_gpu",
        [](const py::array_t<double>& A, const int64_t num_gpus_per_node) -> py::list {
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

#if defined(USE_CUDA)

            gpu::diagonalizeMatrix(evecs, evals, static_cast<int64_t>(nrows_A));

            // TODO: diagonalizeMatrixMultiGPU

#elif defined(USE_HIP)

            if (dim < 8192)
            {
                gpu::diagonalizeMatrix(evecs, evals, static_cast<int64_t>(nrows_A));
            }
            else
            {
                gpu::diagonalizeMatrixMultiGPU(evecs, evals, static_cast<int64_t>(nrows_A), num_gpus_per_node);
            }

#endif

            py::list result;

            result.append(eigenValues);

            result.append(eigenVectors);

            return result;
        },
        "Diagonalizes matrix using GPU.");

    m.def("integrate_vxc_fock_gpu", &gpu::integrateVxcFock, "Integrates Vxc matrix using GPU.");

    m.def("integrate_fxc_fock_gpu", &gpu::integrateFxcFock, "Integrates Fxc matrix using GPU.");

    m.def("integrate_vxc_gradient_gpu", &gpu::integrateVxcGradient, "Integrates Vxc gradient using GPU.");

    m.def(
        "compute_fock_gpu",
        [](const CMolecule&            molecule,
           const CMolecularBasis&      basis,
           const CAODensityMatrix&     densityMatrix,
           const double                prefac_coulomb,
           const std::vector<double>&  frac_exact_exchange_values,
           const std::vector<double>&  omega_values,
           const std::string&          flag_K,
           const double                eri_threshold,
           const double                prelink_threshold,
           const py::array_t<int32_t>& Q_prime_row_indices,
           const py::array_t<int32_t>& Q_prime_col_indices,
           CScreeningData&             screening) -> CDenseMatrix {

            std::string errsize("compute_fock_gpu: Mismatch in the size of Q_prime indices");
            errors::assertMsgCritical(Q_prime_row_indices.size() == Q_prime_col_indices.size(), errsize);

            std::string errstyle("compute_fock_gpu: Expecting contiguous Q_prime indices");
            auto c_style_row = py::detail::check_flags(Q_prime_row_indices.ptr(), py::array::c_style);
            auto c_style_col = py::detail::check_flags(Q_prime_col_indices.ptr(), py::array::c_style);
            errors::assertMsgCritical(c_style_row && c_style_col, errstyle);

            const auto Q_prime_ind_count = static_cast<int64_t>(Q_prime_row_indices.size());
            const auto Q_prime_row_ptr = Q_prime_row_indices.data();
            const auto Q_prime_col_ptr = Q_prime_col_indices.data();

            return gpu::computeFockOnGPU(molecule, basis, densityMatrix, prefac_coulomb,
                                         frac_exact_exchange_values, omega_values, flag_K,
                                         eri_threshold, prelink_threshold, Q_prime_row_ptr,
                                         Q_prime_col_ptr, Q_prime_ind_count, screening);
        },
        "Computes Fock matrix using GPU.");

    m.def("compute_fock_gradient_gpu", &gpu::computeFockGradientOnGPU, "Computes Fock gradient using GPU.");

    m.def("transform_density", &gpu::transformDensity, "Transforms density matrix (spherical to Cartesian).");

    m.def("compute_mixed_basis_overlap_integrals_gpu", &gpu::computeMixedBasisOverlapIntegralsOnGPU, "Computes mixed-basis overlap integral matrix using GPU.");

    m.def("compute_overlap_and_kinetic_energy_integrals_gpu", &gpu::computeOverlapAndKineticEnergyIntegralsOnGPU, "Computes one-electron integral matrices using GPU.");

    m.def("compute_overlap_gradient_gpu", &gpu::computeOverlapGradientOnGPU, "Computes overlap integral gradient using GPU.");

    m.def("compute_kinetic_energy_gradient_gpu", &gpu::computeKineticEnergyGradientOnGPU, "Computes kinetic energy integral gradient using GPU.");

    m.def("compute_nuclear_potential_gradient_gpu", &gpu::computeNuclearPotentialGradientOnGPU, "Computes nuclear potential integral gradient using GPU.");

    m.def("compute_nuclear_potential_integrals_gpu", &gpu::computeNuclearPotentialIntegralsOnGPU, "Computes one-electron integral matrices using GPU.");

    m.def(
        "compute_point_charges_integrals_gpu",
        [](const CMolecule& molecule,
           const CMolecularBasis& basis,
           const CScreeningData& screener,
           const py::array_t<double>& point_charges) -> CDenseMatrix {
            std::string errshape("compute_point_charges_integrals_gpu: Invalid shape of point_charges");
            std::string errstyle("compute_point_charges_integrals_gpu: Expecting contiguous numpy array");
            const auto  ndim    = static_cast<int64_t>(point_charges.shape(0));
            const auto  npoints = static_cast<int64_t>(point_charges.shape(1));
            auto        c_style = py::detail::check_flags(point_charges.ptr(), py::array::c_style);
            errors::assertMsgCritical(ndim == 4, errshape);
            errors::assertMsgCritical(c_style, errstyle);
            return gpu::computePointChargesIntegralsOnGPU(molecule, basis, screener, point_charges.data(), npoints);
        },
        "Computes point charges integrals using GPU.");

    m.def(
        "compute_point_charges_gradient_gpu",
        [](const CMolecule& molecule,
           const CMolecularBasis& basis,
           const CGradientScreeningData& screener,
           const py::array_t<double>& point_charges) -> CDenseMatrix {
            std::string errshape("compute_point_charges_gradient_gpu: Invalid shape of point_charges");
            std::string errstyle("compute_point_charges_gradient_gpu: Expecting contiguous numpy array");
            const auto  ndim    = static_cast<int64_t>(point_charges.shape(0));
            const auto  npoints = static_cast<int64_t>(point_charges.shape(1));
            auto        c_style = py::detail::check_flags(point_charges.ptr(), py::array::c_style);
            errors::assertMsgCritical(ndim == 4, errshape);
            errors::assertMsgCritical(c_style, errstyle);
            return gpu::computePointChargesGradientOnGPU(molecule, basis, screener, point_charges.data(), npoints);
        },
        "Computes point charges integrals using GPU.");

    m.def(
        "compute_electric_dipole_integrals_gpu", &gpu::computeElectricDipoleIntegralsOnGPU, "Computes electric dipole integral matrices using GPU.");

    m.def(
        "compute_linear_momentum_integrals_gpu", &gpu::computeLinearMomentumIntegralsOnGPU, "Computes linear momentum integral matrices using GPU.");

    m.def("compute_angular_momentum_integrals_gpu",
          &gpu::computeAngularMomentumIntegralsOnGPU,
          "Computes angular momentum integral matrices using GPU.");

    m.def("compute_q_matrix_gpu", &gpu::computeQMatrixOnGPU, "Computes Q matrix using GPU.");
}

}  // namespace vlx_gpu
