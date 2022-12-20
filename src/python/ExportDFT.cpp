//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include "ExportDFT.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "NewFunctionalParser.hpp"
#include "XCFuncType.hpp"
#include "XCNewFunctional.hpp"
#include "XCNewIntegrator.hpp"
#include "XCNewMolecularGradient.hpp"
#include "XCPairDensityFunctional.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

CAOKohnShamMatrix
integrate_vxc_pdft(const CXCNewIntegrator&    self,
                   const CAODensityMatrix&    aoDensityMatrix,
                   const py::array_t<double>& Active2DM,
                   const py::array_t<double>& ActiveMOs,
                   const CMolecule&           molecule,
                   const CMolecularBasis&     basis,
                   const CMolecularGrid&      molecularGrid,
                   const std::string&         xcFuncLabel)
{
    // Active2DM

    // check dimension

    std::string errdim("integrate_vxc_pdft, Active2DM: Expecting a 4D numpy array");

    errors::assertMsgCritical(Active2DM.ndim() == 4, errdim);

    // check that the numpy array is c-style contiguous

    std::string errsrc("integrate_vxc_pdft, Active2DM: Expecting a contiguous numpy array in C ordering");

    auto c_style = py::detail::check_flags(Active2DM.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    // Form 4D tensor

    auto nActive = static_cast<int32_t>(Active2DM.shape(0));

    bool same_size =
        ((Active2DM.shape(0) == Active2DM.shape(1)) && (Active2DM.shape(0) == Active2DM.shape(2)) && (Active2DM.shape(0) == Active2DM.shape(3)));

    std::string errsizes("integrate_vxc_pdft, Active2DM: Expecting 4 identical dimensions");

    errors::assertMsgCritical(same_size, errsizes);

    std::vector<double> vec(Active2DM.data(), Active2DM.data() + Active2DM.size());

    CDense4DTensor Tensor_2DM(vec, nActive, nActive, nActive, nActive);

    // active MO

    // Check dimensions

    errdim = "integrate_vxc_pdft, ActiveMOs: Expecting a 2D numpy array";

    errors::assertMsgCritical(ActiveMOs.ndim() == 2, errdim);

    // check that the numpy array is c-style contiguous

    errsrc = "integrate_vxc_pdft, ActiveMOs: Expecting a contiguous numpy array in C ordering";

    c_style = py::detail::check_flags(ActiveMOs.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    auto naos = ActiveMOs.shape(1);

    std::vector<double> vec2(ActiveMOs.data(), ActiveMOs.data() + ActiveMOs.size());

    CDenseMatrix Dense_activeMO(vec2, nActive, naos);

    CAOKohnShamMatrix mat_Vxc(naos, naos, true);

    mat_Vxc.zero();

    CDense4DTensor TwoBodyGradient(naos, nActive, nActive, nActive);

    TwoBodyGradient.zero();

    self.integrateVxcPDFT(mat_Vxc, TwoBodyGradient, molecule, basis, aoDensityMatrix, Tensor_2DM, Dense_activeMO, molecularGrid, xcFuncLabel);

    auto xcene = mat_Vxc.getExchangeCorrelationEnergy();

    return mat_Vxc;
}

// Exports classes/functions in src/dft to python

void
export_dft(py::module& m)
{
    // xcfun enum class

    // clang-format off
    py::enum_<xcfun>(m, "xcfun")
        .value("lda", xcfun::lda)
        .value("gga", xcfun::gga)
        .value("mgga", xcfun::mgga);
    // clang-format on

    // CAOKohnShamMatrix class

    PyClass<CAOKohnShamMatrix>(m, "AOKohnShamMatrix")
        .def(py::init<>())
        .def(py::init<int32_t, int32_t, bool>(), "nrows"_a, "ncols"_a, "is_rest"_a)
        .def("__str__", &CAOKohnShamMatrix::getString)
        .def("get_matrix", &CAOKohnShamMatrix::getReferenceToKohnSham, "Gets constant reference to specific Kohn-Sham matrix.", "beta"_a = false)
        .def(
            "reduce_sum",
            [](CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.reduce_sum(rank, nodes, *comm);
            },
            "Performs reduce_sum for AOKohnShamMatrix object",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a)
        .def(
            "collect",
            [](CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm, int32_t source) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.collect(rank, nodes, *comm, source);
            },
            "Collects AOKohnShamMatrix object.",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a,
            "source"_a)
        .def("get_electrons", &CAOKohnShamMatrix::getNumberOfElectrons, "Gets number of electrons obtained by integrating Kohn-Sham matrix.")
        .def("get_energy", &CAOKohnShamMatrix::getExchangeCorrelationEnergy, "Gets exchange-correlation energy associated with Kohn-Sham matrix.")
        .def(py::self == py::self);

    // XCComponent class
    PyClass<CXCComponent>(m, "XCComponent")
        .def(py::init<const std::string&, const double>(), "label"_a, "coeff"_a)
        .def(py::init<const CXCComponent&>())
        .def("get_scaling_factor", &CXCComponent::getScalingFactor, "Gets scaling factor of XC functional component.")
        .def("get_label", &CXCComponent::getLabel, "Gets name of XC functional component.")
        .def(py::self == py::self);

    // XCNewFunctional class
    PyClass<CXCNewFunctional>(m, "XCNewFunctional")
        .def(py::init<const std::string&, const std::vector<std::string>&, const std::vector<double>&, const double>(),
             "name_of_functional"_a,
             "labels"_a,
             "coeffs"_a,
             "fraction_of_exact_exchange"_a = 0.0)
        .def(py::init<const CXCNewFunctional&>())
        .def(py::self == py::self)
        .def("is_hybrid", &CXCNewFunctional::isHybrid, "Determines whether the XC functional is hybrid.")
        .def("is_undefined", &CXCNewFunctional::isUndefined, "Determines whether the XC function is undefined.")
        .def("get_func_type", &CXCNewFunctional::getFunctionalType, "Gets type of XC functional.")
        .def("get_func_label", &CXCNewFunctional::getFunctionalLabel, "Gets name of XC functional.")
        .def("get_frac_exact_exchange",
             &CXCNewFunctional::getFractionOfExactExchange,
             "Gets fraction of exact Hartree-Fock exchange in XC functional.")
        .def(
            "compute_exc_vxc_for_lda",
            [](const CXCNewFunctional& self, const py::array_t<double>& rho) -> py::list {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_exc_vxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_exc_vxc_for_lda: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                self.compute_exc_vxc_for_lda(npoints, rho.data(), exc.values(), vrho.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfElements()));
                return ret;
            },
            "Computes Exc and Vxc for LDA.",
            "rho"_a)
        .def(
            "compute_exc_vxc_for_gga",
            [](const CXCNewFunctional& self, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::list {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style,
                                          std::string("compute_exc_vxc_for_gga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical((rho_size == npoints * 2) && (sigma_size == npoints * 3),
                                          std::string("compute_exc_vxc_for_gga: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                CDenseMatrix vsigma(npoints, 3);
                self.compute_exc_vxc_for_gga(npoints, rho.data(), sigma.data(), exc.values(), vrho.values(), vsigma.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vsigma.values(), vsigma.getNumberOfElements()));
                return ret;
            },
            "Computes Exc and Vxc for GGA.",
            "rho"_a,
            "sigma"_a);

    // XCPairDensityFunctional class
    PyClass<CXCPairDensityFunctional>(m, "XCPairDensityFunctional")
        .def(py::init<const std::string&, const std::vector<std::string>&, const std::vector<double>&>(),
             "name_of_functional"_a,
             "labels"_a,
             "coeffs"_a)
        .def(py::init<const CXCPairDensityFunctional&>())
        .def(py::self == py::self)
        .def("get_func_label", &CXCPairDensityFunctional::getFunctionalLabel, "Gets functional name.")
        .def("get_func_type", &CXCPairDensityFunctional::getFunctionalType, "Gets functional type.")
        .def(
            "compute_exc_vxc_for_plda",
            [](const CXCPairDensityFunctional& self, const py::array_t<double>& rho) -> py::list {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_exc_vxc_for_plda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_exc_vxc_for_plda: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                self.compute_exc_vxc_for_plda(npoints, rho.data(), exc.values(), vrho.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfElements()));
                return ret;
            },
            "Computes Exc and Vxc for pair-density LDA.",
            "rho"_a);

    // CMolecularGrid class

    PyClass<CMolecularGrid>(m, "MolecularGrid")
        .def(py::init<>())
        .def(py::init<const CMolecularGrid&>())
        .def("number_of_points", &CMolecularGrid::getNumberOfGridPoints)
        .def(
            "x_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesX(), self.getNumberOfGridPoints());
            },
            "Gets X coordinates of grid as numpy array.")
        .def(
            "y_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesY(), self.getNumberOfGridPoints());
            },
            "Gets Y coordinates of grid as numpy array.")
        .def(
            "z_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), self.getNumberOfGridPoints());
            },
            "Gets Z coordinates of grid as numpy array.")
        .def(
            "w_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getWeights(), self.getNumberOfGridPoints());
            },
            "Gets weights of grid as numpy array.")
        .def(
            "broadcast",
            [](CMolecularGrid& self, int32_t rank, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.broadcast(rank, *comm);
            },
            "Broadcasts MolecularGrid object.",
            "rank"_a,
            "py_comm"_a)
        .def("partition_grid_points", &CMolecularGrid::partitionGridPoints)
        .def(
            "distribute_counts_and_displacements",
            [](CMolecularGrid& self, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.distributeCountsAndDisplacements(rank, nodes, *comm);
            },
            "Distributes MolecularGrid counts and displacements.",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a)
        .def(
            "re_distribute_counts_and_displacements",
            [](CMolecularGrid& self, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.reDistributeCountsAndDisplacements(rank, nodes, *comm);
            },
            "Redo distributing MolecularGrid counts and displacements.",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a)
        .def(py::self == py::self);

    // CGridDriver class

    PyClass<CGridDriver>(m, "GridDriver")
        .def(py::init(&vlx_general::create<CGridDriver>), "comm"_a = py::none())
        .def("generate",
             &CGridDriver::generate,
             "Generates molecular grid for molecule. Errors are printed to output stream. Grid generation is distributed within domain of MPI "
             "communicator.",
             "molecule"_a)
        .def("set_level",
             &CGridDriver::setLevel,
             "Sets accuracy level for grid generation. Level: 1-6, where 1 is coarse grid, 5 is ultrafine grid, 6 special benchmarking grid.",
             "gridLevel"_a);

    // CXCNewIntegrator class

    PyClass<CXCNewIntegrator>(m, "XCNewIntegrator")
        .def(py::init(&vlx_general::create<CXCNewIntegrator>), "comm"_a = py::none())
        .def("integrate_vxc_fock",
             &CXCNewIntegrator::integrateVxcFock,
             "Integrates 1st-order exchange-correlation contribution to Kohn-Sham matrix.",
             "molecule"_a,
             "basis"_a,
             "densityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate_fxc_fock",
             &CXCNewIntegrator::integrateFxcFock,
             "Integrates 2nd-order exchange-correlation contribution to Fock matrix.",
             "aoFockMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "rwDensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate_kxc_fock",
             &CXCNewIntegrator::integrateKxcFock,
             "Integrates 3rd-order exchange-correlation contribution to Fock matrix.",
             "aoFockMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "rwDensityMatrix"_a,
             "rw2DensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a,
             "quadMode"_a)
        .def("integrate_vxc_pdft", &integrate_vxc_pdft)
        .def(
            "compute_gto_values",
            [](CXCNewIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid)
                -> py::array_t<double> {
                auto gtovalues = self.computeGtoValuesOnGridPoints(molecule, basis, molecularGrid);
                return vlx_general::pointer_to_numpy(gtovalues.values(), gtovalues.getNumberOfRows(), gtovalues.getNumberOfColumns());
            },
            "Computes GTO values on grid points.",
            "molecule"_a,
            "basis"_a,
            "molecularGrid"_a)
        .def(
            "compute_gto_values_and_derivatives",
            [](CXCNewIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> py::list {
                auto     gtovaluesderivs = self.computeGtoValuesAndDerivativesOnGridPoints(molecule, basis, molecularGrid);
                py::list ret;
                for (int32_t i = 0; i < static_cast<int32_t>(gtovaluesderivs.size()); i++)
                {
                    ret.append(vlx_general::pointer_to_numpy(
                        gtovaluesderivs[i].values(), gtovaluesderivs[i].getNumberOfRows(), gtovaluesderivs[i].getNumberOfColumns()));
                }
                return ret;
            },
            "Computes GTO values and derivatives on grid points.",
            "molecule"_a,
            "basis"_a,
            "molecularGrid"_a)
        .def(
            "compute_gto_values_and_derivatives",
            [](CXCNewIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const py::array_t<double>& points) -> py::list {
                auto points_c_style = py::detail::check_flags(points.ptr(), py::array::c_style);
                errors::assertMsgCritical(points_c_style,
                                          std::string("compute_gto_values_and_derivatives_on_points: Expecting C-style contiguous numpy array"));
                errors::assertMsgCritical(points.shape(0) == 3,
                                          std::string("compute_gto_values_and_derivatives_on_points: Expecting numpy array of shape (3,N)"));
                auto     npoints         = static_cast<int32_t>(points.shape(1));
                auto     xcoords         = points.data();
                auto     ycoords         = xcoords + npoints;
                auto     zcoords         = ycoords + npoints;
                auto     gtovaluesderivs = self.computeGtoValuesAndDerivativesOnGridPoints(molecule, basis, npoints, xcoords, ycoords, zcoords);
                py::list ret;
                for (int32_t i = 0; i < static_cast<int32_t>(gtovaluesderivs.size()); i++)
                {
                    ret.append(vlx_general::pointer_to_numpy(
                        gtovaluesderivs[i].values(), gtovaluesderivs[i].getNumberOfRows(), gtovaluesderivs[i].getNumberOfColumns()));
                }
                return ret;
            },
            "Computes GTO values and derivatives on grid points.",
            "molecule"_a,
            "basis"_a,
            "molecularGrid"_a)
        .def(
            "compute_exc_vxc_for_lda",
            [](CXCNewIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::list {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_exc_vxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_exc_vxc_for_lda: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                self.computeExcVxcForLDA(xcFuncLabel, npoints, rho.data(), exc.values(), vrho.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfElements()));
                return ret;
            },
            "Computes Exc and Vxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_exc_vxc_for_gga",
            [](CXCNewIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::list {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style,
                                          std::string("compute_exc_vxc_for_gga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical((rho_size == npoints * 2) && (sigma_size == npoints * 3),
                                          std::string("compute_exc_vxc_for_gga: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                CDenseMatrix vsigma(npoints, 3);
                self.computeExcVxcForGGA(xcFuncLabel, npoints, rho.data(), sigma.data(), exc.values(), vrho.values(), vsigma.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(vsigma.values(), vsigma.getNumberOfElements()));
                return ret;
            },
            "Computes Exc and Vxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a)
        .def(
            "compute_fxc_for_lda",
            [](CXCNewIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::array_t<double> {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_fxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_fxc_for_lda: Inconsistent array size"));
                CDenseMatrix v2rho2(npoints, 3);
                self.computeFxcForLDA(xcFuncLabel, npoints, rho.data(), v2rho2.values());
                return vlx_general::pointer_to_numpy(v2rho2.values(), v2rho2.getNumberOfElements());
            },
            "Computes Fxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_fxc_for_gga",
            [](CXCNewIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::list {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style, std::string("compute_fxc_for_gga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical((rho_size == npoints * 2) && (sigma_size == npoints * 3),
                                          std::string("compute_fxc_for_gga: Inconsistent array size"));
                CDenseMatrix v2rho2(npoints, 3);
                CDenseMatrix v2rhosigma(npoints, 6);
                CDenseMatrix v2sigma2(npoints, 6);
                self.computeFxcForGGA(xcFuncLabel, npoints, rho.data(), sigma.data(), v2rho2.values(), v2rhosigma.values(), v2sigma2.values());
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(v2rho2.values(), v2rho2.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(v2rhosigma.values(), v2rhosigma.getNumberOfElements()));
                ret.append(vlx_general::pointer_to_numpy(v2sigma2.values(), v2sigma2.getNumberOfElements()));
                return ret;
            },
            "Computes Fxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a);

    // CXCNewMolecularGradient class

    PyClass<CXCNewMolecularGradient>(m, "XCNewMolecularGradient")
        .def(py::init(&vlx_general::create<CXCNewMolecularGradient>), "comm"_a = py::none())
        .def(
            "integrate_vxc_gradient",
            [](CXCNewMolecularGradient& self,
               const CMolecule&         molecule,
               const CMolecularBasis&   basis,
               const CAODensityMatrix&  gsDensityMatrix,
               const CMolecularGrid&    molecularGrid,
               const std::string&       xcFuncLabel) -> py::array_t<double> {
                auto molgrad = self.integrateVxcGradient(molecule, basis, gsDensityMatrix, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), molgrad.getNumberOfRows(), molgrad.getNumberOfColumns());
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_vxc_gradient",
            [](CXCNewMolecularGradient& self,
               const CMolecule&         molecule,
               const CMolecularBasis&   basis,
               const CAODensityMatrix&  rwDensityMatrix,
               const CAODensityMatrix&  gsDensityMatrix,
               const CMolecularGrid&    molecularGrid,
               const std::string&       xcFuncLabel) -> py::array_t<double> {
                auto molgrad = self.integrateVxcGradient(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), molgrad.getNumberOfRows(), molgrad.getNumberOfColumns());
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityMatrix"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_fxc_gradient",
            [](CXCNewMolecularGradient& self,
               const CMolecule&         molecule,
               const CMolecularBasis&   basis,
               const CAODensityMatrix&  rwDensityMatrixOne,
               const CAODensityMatrix&  rwDensityMatrixTwo,
               const CAODensityMatrix&  gsDensityMatrix,
               const CMolecularGrid&    molecularGrid,
               const std::string&       xcFuncLabel) -> py::array_t<double> {
                auto molgrad =
                    self.integrateFxcGradient(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo, gsDensityMatrix, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), molgrad.getNumberOfRows(), molgrad.getNumberOfColumns());
            },
            "Integrates 2nd-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityMatrixOne"_a,
            "rwDensityMatrixTwo"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_kxc_gradient",
            [](CXCNewMolecularGradient& self,
               const CMolecule&         molecule,
               const CMolecularBasis&   basis,
               const CAODensityMatrix&  rwDensityMatrixOne,
               const CAODensityMatrix&  rwDensityMatrixTwo,
               const CAODensityMatrix&  gsDensityMatrix,
               const CMolecularGrid&    molecularGrid,
               const std::string&       xcFuncLabel) -> py::array_t<double> {
                auto molgrad =
                    self.integrateKxcGradient(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo, gsDensityMatrix, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), molgrad.getNumberOfRows(), molgrad.getNumberOfColumns());
            },
            "Integrates 3rd-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityMatrixOne"_a,
            "rwDensityMatrixTwo"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a);

    // CDensityGrid class

    PyClass<CDensityGrid>(m, "DensityGrid")
        .def(py::init<>())
        .def(py::init<const CDensityGrid&>())
        .def("number_of_points", &CDensityGrid::getNumberOfGridPoints)
        .def("number_of_density_matrices", &CDensityGrid::getNumberOfDensityMatrices)
        .def(
            "density_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensity(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "density_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensity(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_norm_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient norm on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_norm_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient norm on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_product_ab_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.mixedDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets mixed density gradient product on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_x_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientX(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient X component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_y_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientY(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient Y component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_z_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientZ(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient Z component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_x_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientX(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient X component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_y_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientY(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient Y component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_z_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientZ(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient Z component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(py::self == py::self);

    // exposing functions

    m.def("to_xcfun", &to_xcfun, "Converts string label to its enumerate class value.", "label"_a);

    m.def("available_functionals", &newvxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");

    m.def("new_parse_xc_func",
          &newvxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
