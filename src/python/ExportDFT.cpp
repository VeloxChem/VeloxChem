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
#include "Dense4DTensor.hpp"
#include "DensityGrid.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "FunctionalParser.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"
#include "XCFunctional.hpp"
#include "XCIntegrator.hpp"
#include "XCMolecularGradient.hpp"
#include "XCMolecularHessian.hpp"
#include "XCPairDensityFunctional.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

py::list
integrate_vxc_pdft(const CXCIntegrator&       self,
                   const CAODensityMatrix&    aoDensityMatrix,
                   const py::array_t<double>& active2DM,
                   const py::array_t<double>& activeMOs,
                   const CMolecule&           molecule,
                   const CMolecularBasis&     basis,
                   const CMolecularGrid&      molecularGrid,
                   const std::string&         xcFuncLabel)
{
    // active2DM

    // check dimension

    std::string errdim("integrate_vxc_pdft, active2DM: Expecting a 4D numpy array");

    errors::assertMsgCritical(active2DM.ndim() == 4, errdim);

    // check that the numpy array is c-style contiguous

    std::string errsrc("integrate_vxc_pdft, active2DM: Expecting a C-style contiguous numpy array");

    auto c_style = py::detail::check_flags(active2DM.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    // Form 4D tensor

    auto n_active = static_cast<int32_t>(active2DM.shape(0));

    bool same_size =
        ((active2DM.shape(0) == active2DM.shape(1)) && (active2DM.shape(0) == active2DM.shape(2)) && (active2DM.shape(0) == active2DM.shape(3)));

    std::string errsizes("integrate_vxc_pdft, active2DM: Expecting 4 identical dimensions");

    errors::assertMsgCritical(same_size, errsizes);

    std::vector<double> vec(active2DM.data(), active2DM.data() + active2DM.size());

    CDense4DTensor tensor2DM(vec, n_active, n_active, n_active, n_active);

    // activeMOs

    // Check dimensions

    errdim = "integrate_vxc_pdft, activeMOs: Expecting a 2D numpy array";

    errors::assertMsgCritical(activeMOs.ndim() == 2, errdim);

    // check that the numpy array is c-style contiguous

    errsrc = "integrate_vxc_pdft, activeMOs: Expecting a C-style contiguous numpy array";

    c_style = py::detail::check_flags(activeMOs.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    auto naos = activeMOs.shape(1);

    std::vector<double> vec2(activeMOs.data(), activeMOs.data() + activeMOs.size());

    CDenseMatrix denseActiveMO(vec2, n_active, naos);

    // Create output tensors

    CAOKohnShamMatrix matrixVxc(naos, naos, true);

    matrixVxc.zero();

    CDense4DTensor tensorWxc(naos, n_active, n_active, n_active);

    tensorWxc.zero();

    self.integrateVxcPDFT(matrixVxc, tensorWxc, molecule, basis, aoDensityMatrix, tensor2DM, denseActiveMO, molecularGrid, xcFuncLabel);

    py::list returnList;

    returnList.append(matrixVxc);

    returnList.append(vlx_general::pointer_to_numpy(tensorWxc.values(), naos, n_active * n_active * n_active));

    return returnList;
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

    // XCFunctional class
    PyClass<CXCFunctional>(m, "XCFunctional")
        .def(py::init<const std::string&, const std::vector<std::string>&, const std::vector<double>&, const double>(),
             "name_of_functional"_a,
             "labels"_a,
             "coeffs"_a,
             "fraction_of_exact_exchange"_a = 0.0)
        .def(py::init<const CXCFunctional&>())
        .def(py::self == py::self)
        .def("is_hybrid", &CXCFunctional::isHybrid, "Determines whether the XC functional is hybrid.")
        .def("is_undefined", &CXCFunctional::isUndefined, "Determines whether the XC function is undefined.")
        .def("get_func_type", &CXCFunctional::getFunctionalType, "Gets type of XC functional.")
        .def("get_func_label", &CXCFunctional::getFunctionalLabel, "Gets name of XC functional.")
        .def("get_frac_exact_exchange", &CXCFunctional::getFractionOfExactExchange, "Gets fraction of exact Hartree-Fock exchange in XC functional.")
        .def(
            "compute_exc_vxc_for_lda",
            [](const CXCFunctional& self, const py::array_t<double>& rho) -> py::list {
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
            [](const CXCFunctional& self, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::list {
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

    // CXCIntegrator class

    PyClass<CXCIntegrator>(m, "XCIntegrator")
        .def(py::init(&vlx_general::create<CXCIntegrator>), "comm"_a = py::none())
        .def("integrate_vxc_fock",
             &CXCIntegrator::integrateVxcFock,
             "Integrates 1st-order exchange-correlation contribution to Kohn-Sham matrix.",
             "molecule"_a,
             "basis"_a,
             "densityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate_fxc_fock",
             &CXCIntegrator::integrateFxcFock,
             "Integrates 2nd-order exchange-correlation contribution to Fock matrix.",
             "aoFockMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "rwDensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate_kxc_fock",
             &CXCIntegrator::integrateKxcFock,
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
        .def("integrate_lxc_fock",
             &CXCIntegrator::integrateLxcFock,
             "Integrates 4rd-order exchange-correlation contribution to Fock matrix.",
             "aoFockMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "rwDensityMatrix"_a,
             "rw2DensityMatrix"_a,
             "rw3DensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a,
             "cubeMode"_a)
        .def("integrate_kxclxc_fock",
             &CXCIntegrator::integrateKxcLxcFock,
             "Integrates 4rd-order exchange-correlation contribution to Fock matrix.",
             "aoFockMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "rwDensityMatrix"_a,
             "rw2DensityMatrix"_a,
             "rw3DensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a,
             "cubeMode"_a)
        .def("integrate_vxc_pdft", &integrate_vxc_pdft)
        .def(
            "compute_gto_values",
            [](CXCIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid)
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
            [](CXCIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> py::list {
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
            [](CXCIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const py::array_t<double>& points) -> py::list {
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
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::dict {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_exc_vxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_exc_vxc_for_lda: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_exc_vxc_for_lda(npoints, rho.data(), exc.values(), vrho.values());
                py::dict ret;
                ret["exc"]  = vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfRows(), exc.getNumberOfColumns());
                ret["vrho"] = vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfRows(), vrho.getNumberOfColumns());
                return ret;
            },
            "Computes Exc and Vxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_fxc_for_lda",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::dict {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_fxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_fxc_for_lda: Inconsistent array size"));
                CDenseMatrix v2rho2(npoints, 3);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_fxc_for_lda(npoints, rho.data(), v2rho2.values());
                py::dict ret;
                ret["v2rho2"] = vlx_general::pointer_to_numpy(v2rho2.values(), v2rho2.getNumberOfRows(), v2rho2.getNumberOfColumns());
                return ret;
            },
            "Computes Fxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_kxc_for_lda",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::dict {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_kxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_kxc_for_lda: Inconsistent array size"));
                CDenseMatrix v3rho3(npoints, 4);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_kxc_for_lda(npoints, rho.data(), v3rho3.values());
                py::dict ret;
                ret["v3rho3"] = vlx_general::pointer_to_numpy(v3rho3.values(), v3rho3.getNumberOfRows(), v3rho3.getNumberOfColumns());
                return ret;
            },
            "Computes Kxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_lxc_for_lda",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho) -> py::dict {
                auto rho_c_style = py::detail::check_flags(rho.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style, std::string("compute_lxc_for_lda: Expecting C-style contiguous numpy array"));
                auto rho_size = static_cast<int32_t>(rho.size());
                auto npoints  = rho_size / 2;
                errors::assertMsgCritical(rho_size == npoints * 2, std::string("compute_lxc_for_lda: Inconsistent array size"));
                CDenseMatrix v4rho4(npoints, 5);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_lxc_for_lda(npoints, rho.data(), v4rho4.values());
                py::dict ret;
                ret["v4rho4"] = vlx_general::pointer_to_numpy(v4rho4.values(), v4rho4.getNumberOfRows(), v4rho4.getNumberOfColumns());
                return ret;
            },
            "Computes Lxc for LDA.",
            "xcFuncLabel"_a,
            "rho"_a)
        .def(
            "compute_exc_vxc_for_gga",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::dict {
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
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_exc_vxc_for_gga(npoints, rho.data(), sigma.data(), exc.values(), vrho.values(), vsigma.values());
                py::dict ret;
                ret["exc"]    = vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfRows(), exc.getNumberOfColumns());
                ret["vrho"]   = vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfRows(), vrho.getNumberOfColumns());
                ret["vsigma"] = vlx_general::pointer_to_numpy(vsigma.values(), vsigma.getNumberOfRows(), vsigma.getNumberOfColumns());
                return ret;
            },
            "Computes Exc and Vxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a)
        .def(
            "compute_fxc_for_gga",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::dict {
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
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_fxc_for_gga(npoints, rho.data(), sigma.data(), v2rho2.values(), v2rhosigma.values(), v2sigma2.values());
                py::dict ret;
                ret["v2rho2"]     = vlx_general::pointer_to_numpy(v2rho2.values(), v2rho2.getNumberOfRows(), v2rho2.getNumberOfColumns());
                ret["v2rhosigma"] = vlx_general::pointer_to_numpy(v2rhosigma.values(), v2rhosigma.getNumberOfRows(), v2rhosigma.getNumberOfColumns());
                ret["v2sigma2"]   = vlx_general::pointer_to_numpy(v2sigma2.values(), v2sigma2.getNumberOfRows(), v2sigma2.getNumberOfColumns());
                return ret;
            },
            "Computes Fxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a)
        .def(
            "compute_kxc_for_gga",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style, std::string("compute_kxc_for_gga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical((rho_size == npoints * 2) && (sigma_size == npoints * 3),
                                          std::string("compute_kxc_for_gga: Inconsistent array size"));
                CDenseMatrix v3rho3(npoints, 4);
                CDenseMatrix v3rho2sigma(npoints, 9);
                CDenseMatrix v3rhosigma2(npoints, 12);
                CDenseMatrix v3sigma3(npoints, 10);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_kxc_for_gga(
                    npoints, rho.data(), sigma.data(), v3rho3.values(), v3rho2sigma.values(), v3rhosigma2.values(), v3sigma3.values());
                py::dict ret;
                ret["v3rho3"] = vlx_general::pointer_to_numpy(v3rho3.values(), v3rho3.getNumberOfRows(), v3rho3.getNumberOfColumns());
                ret["v3rho2sigma"] =
                    vlx_general::pointer_to_numpy(v3rho2sigma.values(), v3rho2sigma.getNumberOfRows(), v3rho2sigma.getNumberOfColumns());
                ret["v3rhosigma2"] =
                    vlx_general::pointer_to_numpy(v3rhosigma2.values(), v3rhosigma2.getNumberOfRows(), v3rhosigma2.getNumberOfColumns());
                ret["v3sigma3"] = vlx_general::pointer_to_numpy(v3sigma3.values(), v3sigma3.getNumberOfRows(), v3sigma3.getNumberOfColumns());
                return ret;
            },
            "Computes Kxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a)
        .def(
            "compute_lxc_for_gga",
            [](CXCIntegrator& self, const std::string& xcFuncLabel, const py::array_t<double>& rho, const py::array_t<double>& sigma) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style, std::string("compute_lxc_for_gga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical((rho_size == npoints * 2) && (sigma_size == npoints * 3),
                                          std::string("compute_lxc_for_gga: Inconsistent array size"));
                CDenseMatrix v4rho4(npoints, 5);
                CDenseMatrix v4rho3sigma(npoints, 12);
                CDenseMatrix v4rho2sigma2(npoints, 18);
                CDenseMatrix v4rhosigma3(npoints, 20);
                CDenseMatrix v4sigma4(npoints, 15);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_lxc_for_gga(npoints,
                                         rho.data(),
                                         sigma.data(),
                                         v4rho4.values(),
                                         v4rho3sigma.values(),
                                         v4rho2sigma2.values(),
                                         v4rhosigma3.values(),
                                         v4sigma4.values());
                py::dict ret;
                ret["v4rho4"] = vlx_general::pointer_to_numpy(v4rho4.values(), v4rho4.getNumberOfRows(), v4rho4.getNumberOfColumns());
                ret["v4rho3sigma"] =
                    vlx_general::pointer_to_numpy(v4rho3sigma.values(), v4rho3sigma.getNumberOfRows(), v4rho3sigma.getNumberOfColumns());
                ret["v4rho2sigma2"] =
                    vlx_general::pointer_to_numpy(v4rho2sigma2.values(), v4rho2sigma2.getNumberOfRows(), v4rho2sigma2.getNumberOfColumns());
                ret["v4rhosigma3"] =
                    vlx_general::pointer_to_numpy(v4rhosigma3.values(), v4rhosigma3.getNumberOfRows(), v4rhosigma3.getNumberOfColumns());
                ret["v4sigma4"] = vlx_general::pointer_to_numpy(v4sigma4.values(), v4sigma4.getNumberOfRows(), v4sigma4.getNumberOfColumns());
                return ret;
            },
            "Computes Lxc for GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a)
        .def(
            "compute_exc_vxc_for_mgga",
            [](CXCIntegrator&             self,
               const std::string&         xcFuncLabel,
               const py::array_t<double>& rho,
               const py::array_t<double>& sigma,
               const py::array_t<double>& lapl,
               const py::array_t<double>& tau) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                auto lapl_c_style  = py::detail::check_flags(lapl.ptr(), py::array::c_style);
                auto tau_c_style   = py::detail::check_flags(tau.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style && lapl_c_style && tau_c_style,
                                          std::string("compute_exc_vxc_for_mgga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto lapl_size  = static_cast<int32_t>(lapl.size());
                auto tau_size   = static_cast<int32_t>(tau.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical(
                    (rho_size == npoints * 2) && (sigma_size == npoints * 3) && (lapl_size == npoints * 2) && (tau_size == npoints * 2),
                    std::string("compute_exc_vxc_for_mgga: Inconsistent array size"));
                CDenseMatrix exc(npoints, 1);
                CDenseMatrix vrho(npoints, 2);
                CDenseMatrix vsigma(npoints, 3);
                CDenseMatrix vlapl(npoints, 2);
                CDenseMatrix vtau(npoints, 2);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_exc_vxc_for_mgga(npoints,
                                              rho.data(),
                                              sigma.data(),
                                              lapl.data(),
                                              tau.data(),
                                              exc.values(),
                                              vrho.values(),
                                              vsigma.values(),
                                              vlapl.values(),
                                              vtau.values());
                py::dict ret;
                ret["exc"]    = vlx_general::pointer_to_numpy(exc.values(), exc.getNumberOfRows(), exc.getNumberOfColumns());
                ret["vrho"]   = vlx_general::pointer_to_numpy(vrho.values(), vrho.getNumberOfRows(), vrho.getNumberOfColumns());
                ret["vsigma"] = vlx_general::pointer_to_numpy(vsigma.values(), vsigma.getNumberOfRows(), vsigma.getNumberOfColumns());
                ret["vlapl"]  = vlx_general::pointer_to_numpy(vlapl.values(), vlapl.getNumberOfRows(), vlapl.getNumberOfColumns());
                ret["vtau"]   = vlx_general::pointer_to_numpy(vtau.values(), vtau.getNumberOfRows(), vtau.getNumberOfColumns());
                return ret;
            },
            "Computes Exc and Vxc for meta-GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a,
            "lapl"_a,
            "tau"_a)
        .def(
            "compute_fxc_for_mgga",
            [](CXCIntegrator&             self,
               const std::string&         xcFuncLabel,
               const py::array_t<double>& rho,
               const py::array_t<double>& sigma,
               const py::array_t<double>& lapl,
               const py::array_t<double>& tau) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                auto lapl_c_style  = py::detail::check_flags(lapl.ptr(), py::array::c_style);
                auto tau_c_style   = py::detail::check_flags(tau.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style && lapl_c_style && tau_c_style,
                                          std::string("compute_fxc_for_mgga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto lapl_size  = static_cast<int32_t>(lapl.size());
                auto tau_size   = static_cast<int32_t>(tau.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical(
                    (rho_size == npoints * 2) && (sigma_size == npoints * 3) && (lapl_size == npoints * 2) && (tau_size == npoints * 2),
                    std::string("compute_fxc_for_mgga: Inconsistent array size"));
                CDenseMatrix v2rho2(npoints, 3);
                CDenseMatrix v2rhosigma(npoints, 6);
                CDenseMatrix v2rholapl(npoints, 4);
                CDenseMatrix v2rhotau(npoints, 4);
                CDenseMatrix v2sigma2(npoints, 6);
                CDenseMatrix v2sigmalapl(npoints, 6);
                CDenseMatrix v2sigmatau(npoints, 6);
                CDenseMatrix v2lapl2(npoints, 3);
                CDenseMatrix v2lapltau(npoints, 4);
                CDenseMatrix v2tau2(npoints, 3);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_fxc_for_mgga(npoints,
                                          rho.data(),
                                          sigma.data(),
                                          lapl.data(),
                                          tau.data(),
                                          v2rho2.values(),
                                          v2rhosigma.values(),
                                          v2rholapl.values(),
                                          v2rhotau.values(),
                                          v2sigma2.values(),
                                          v2sigmalapl.values(),
                                          v2sigmatau.values(),
                                          v2lapl2.values(),
                                          v2lapltau.values(),
                                          v2tau2.values());
                py::dict ret;
                ret["v2rho2"]     = vlx_general::pointer_to_numpy(v2rho2.values(), v2rho2.getNumberOfRows(), v2rho2.getNumberOfColumns());
                ret["v2rhosigma"] = vlx_general::pointer_to_numpy(v2rhosigma.values(), v2rhosigma.getNumberOfRows(), v2rhosigma.getNumberOfColumns());
                ret["v2rholapl"]  = vlx_general::pointer_to_numpy(v2rholapl.values(), v2rholapl.getNumberOfRows(), v2rholapl.getNumberOfColumns());
                ret["v2rhotau"]   = vlx_general::pointer_to_numpy(v2rhotau.values(), v2rhotau.getNumberOfRows(), v2rhotau.getNumberOfColumns());
                ret["v2sigma2"]   = vlx_general::pointer_to_numpy(v2sigma2.values(), v2sigma2.getNumberOfRows(), v2sigma2.getNumberOfColumns());
                ret["v2sigmalapl"] =
                    vlx_general::pointer_to_numpy(v2sigmalapl.values(), v2sigmalapl.getNumberOfRows(), v2sigmalapl.getNumberOfColumns());
                ret["v2sigmatau"] = vlx_general::pointer_to_numpy(v2sigmatau.values(), v2sigmatau.getNumberOfRows(), v2sigmatau.getNumberOfColumns());
                ret["v2lapl2"]    = vlx_general::pointer_to_numpy(v2lapl2.values(), v2lapl2.getNumberOfRows(), v2lapl2.getNumberOfColumns());
                ret["v2lapltau"]  = vlx_general::pointer_to_numpy(v2lapltau.values(), v2lapltau.getNumberOfRows(), v2lapltau.getNumberOfColumns());
                ret["v2tau2"]     = vlx_general::pointer_to_numpy(v2tau2.values(), v2tau2.getNumberOfRows(), v2tau2.getNumberOfColumns());
                return ret;
            },
            "Computes Fxc for meta-GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a,
            "lapl"_a,
            "tau"_a)
        .def(
            "compute_kxc_for_mgga",
            [](CXCIntegrator&             self,
               const std::string&         xcFuncLabel,
               const py::array_t<double>& rho,
               const py::array_t<double>& sigma,
               const py::array_t<double>& lapl,
               const py::array_t<double>& tau) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                auto lapl_c_style  = py::detail::check_flags(lapl.ptr(), py::array::c_style);
                auto tau_c_style   = py::detail::check_flags(tau.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style && lapl_c_style && tau_c_style,
                                          std::string("compute_kxc_for_mgga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto lapl_size  = static_cast<int32_t>(lapl.size());
                auto tau_size   = static_cast<int32_t>(tau.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical(
                    (rho_size == npoints * 2) && (sigma_size == npoints * 3) && (lapl_size == npoints * 2) && (tau_size == npoints * 2),
                    std::string("compute_kxc_for_mgga: Inconsistent array size"));
                CDenseMatrix v3rho3(npoints, 4);
                CDenseMatrix v3rho2sigma(npoints, 9);
                CDenseMatrix v3rho2lapl(npoints, 6);
                CDenseMatrix v3rho2tau(npoints, 6);
                CDenseMatrix v3rhosigma2(npoints, 12);
                CDenseMatrix v3rhosigmalapl(npoints, 12);
                CDenseMatrix v3rhosigmatau(npoints, 12);
                CDenseMatrix v3rholapl2(npoints, 6);
                CDenseMatrix v3rholapltau(npoints, 8);
                CDenseMatrix v3rhotau2(npoints, 6);
                CDenseMatrix v3sigma3(npoints, 10);
                CDenseMatrix v3sigma2lapl(npoints, 12);
                CDenseMatrix v3sigma2tau(npoints, 12);
                CDenseMatrix v3sigmalapl2(npoints, 9);
                CDenseMatrix v3sigmalapltau(npoints, 12);
                CDenseMatrix v3sigmatau2(npoints, 9);
                CDenseMatrix v3lapl3(npoints, 4);
                CDenseMatrix v3lapl2tau(npoints, 6);
                CDenseMatrix v3lapltau2(npoints, 6);
                CDenseMatrix v3tau3(npoints, 4);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_kxc_for_mgga(npoints,
                                          rho.data(),
                                          sigma.data(),
                                          lapl.data(),
                                          tau.data(),
                                          v3rho3.values(),
                                          v3rho2sigma.values(),
                                          v3rho2lapl.values(),
                                          v3rho2tau.values(),
                                          v3rhosigma2.values(),
                                          v3rhosigmalapl.values(),
                                          v3rhosigmatau.values(),
                                          v3rholapl2.values(),
                                          v3rholapltau.values(),
                                          v3rhotau2.values(),
                                          v3sigma3.values(),
                                          v3sigma2lapl.values(),
                                          v3sigma2tau.values(),
                                          v3sigmalapl2.values(),
                                          v3sigmalapltau.values(),
                                          v3sigmatau2.values(),
                                          v3lapl3.values(),
                                          v3lapl2tau.values(),
                                          v3lapltau2.values(),
                                          v3tau3.values());
                py::dict ret;
                ret["v3rho3"] = vlx_general::pointer_to_numpy(v3rho3.values(), v3rho3.getNumberOfRows(), v3rho3.getNumberOfColumns());
                ret["v3rho2sigma"] =
                    vlx_general::pointer_to_numpy(v3rho2sigma.values(), v3rho2sigma.getNumberOfRows(), v3rho2sigma.getNumberOfColumns());
                ret["v3rho2lapl"] = vlx_general::pointer_to_numpy(v3rho2lapl.values(), v3rho2lapl.getNumberOfRows(), v3rho2lapl.getNumberOfColumns());
                ret["v3rho2tau"]  = vlx_general::pointer_to_numpy(v3rho2tau.values(), v3rho2tau.getNumberOfRows(), v3rho2tau.getNumberOfColumns());
                ret["v3rhosigma2"] =
                    vlx_general::pointer_to_numpy(v3rhosigma2.values(), v3rhosigma2.getNumberOfRows(), v3rhosigma2.getNumberOfColumns());
                ret["v3rhosigmalapl"] =
                    vlx_general::pointer_to_numpy(v3rhosigmalapl.values(), v3rhosigmalapl.getNumberOfRows(), v3rhosigmalapl.getNumberOfColumns());
                ret["v3rhosigmatau"] =
                    vlx_general::pointer_to_numpy(v3rhosigmatau.values(), v3rhosigmatau.getNumberOfRows(), v3rhosigmatau.getNumberOfColumns());
                ret["v3rholapl2"] = vlx_general::pointer_to_numpy(v3rholapl2.values(), v3rholapl2.getNumberOfRows(), v3rholapl2.getNumberOfColumns());
                ret["v3rholapltau"] =
                    vlx_general::pointer_to_numpy(v3rholapltau.values(), v3rholapltau.getNumberOfRows(), v3rholapltau.getNumberOfColumns());
                ret["v3rhotau2"] = vlx_general::pointer_to_numpy(v3rhotau2.values(), v3rhotau2.getNumberOfRows(), v3rhotau2.getNumberOfColumns());
                ret["v3sigma3"]  = vlx_general::pointer_to_numpy(v3sigma3.values(), v3sigma3.getNumberOfRows(), v3sigma3.getNumberOfColumns());
                ret["v3sigma2lapl"] =
                    vlx_general::pointer_to_numpy(v3sigma2lapl.values(), v3sigma2lapl.getNumberOfRows(), v3sigma2lapl.getNumberOfColumns());
                ret["v3sigma2tau"] =
                    vlx_general::pointer_to_numpy(v3sigma2tau.values(), v3sigma2tau.getNumberOfRows(), v3sigma2tau.getNumberOfColumns());
                ret["v3sigmalapl2"] =
                    vlx_general::pointer_to_numpy(v3sigmalapl2.values(), v3sigmalapl2.getNumberOfRows(), v3sigmalapl2.getNumberOfColumns());
                ret["v3sigmalapltau"] =
                    vlx_general::pointer_to_numpy(v3sigmalapltau.values(), v3sigmalapltau.getNumberOfRows(), v3sigmalapltau.getNumberOfColumns());
                ret["v3sigmatau2"] =
                    vlx_general::pointer_to_numpy(v3sigmatau2.values(), v3sigmatau2.getNumberOfRows(), v3sigmatau2.getNumberOfColumns());
                ret["v3lapl3"]    = vlx_general::pointer_to_numpy(v3lapl3.values(), v3lapl3.getNumberOfRows(), v3lapl3.getNumberOfColumns());
                ret["v3lapl2tau"] = vlx_general::pointer_to_numpy(v3lapl2tau.values(), v3lapl2tau.getNumberOfRows(), v3lapl2tau.getNumberOfColumns());
                ret["v3lapltau2"] = vlx_general::pointer_to_numpy(v3lapltau2.values(), v3lapltau2.getNumberOfRows(), v3lapltau2.getNumberOfColumns());
                ret["v3tau3"]     = vlx_general::pointer_to_numpy(v3tau3.values(), v3tau3.getNumberOfRows(), v3tau3.getNumberOfColumns());
                return ret;
            },
            "Computes Kxc for meta-GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a,
            "lapl"_a,
            "tau"_a)
        .def(
            "compute_lxc_for_mgga",
            [](CXCIntegrator&             self,
               const std::string&         xcFuncLabel,
               const py::array_t<double>& rho,
               const py::array_t<double>& sigma,
               const py::array_t<double>& lapl,
               const py::array_t<double>& tau) -> py::dict {
                auto rho_c_style   = py::detail::check_flags(rho.ptr(), py::array::c_style);
                auto sigma_c_style = py::detail::check_flags(sigma.ptr(), py::array::c_style);
                auto lapl_c_style  = py::detail::check_flags(lapl.ptr(), py::array::c_style);
                auto tau_c_style   = py::detail::check_flags(tau.ptr(), py::array::c_style);
                errors::assertMsgCritical(rho_c_style && sigma_c_style && lapl_c_style && tau_c_style,
                                          std::string("compute_lxc_for_mgga: Expecting C-style contiguous numpy array"));
                auto rho_size   = static_cast<int32_t>(rho.size());
                auto sigma_size = static_cast<int32_t>(sigma.size());
                auto lapl_size  = static_cast<int32_t>(lapl.size());
                auto tau_size   = static_cast<int32_t>(tau.size());
                auto npoints    = rho_size / 2;
                errors::assertMsgCritical(
                    (rho_size == npoints * 2) && (sigma_size == npoints * 3) && (lapl_size == npoints * 2) && (tau_size == npoints * 2),
                    std::string("compute_lxc_for_mgga: Inconsistent array size"));
                CDenseMatrix v4rho4(npoints, 5);
                CDenseMatrix v4rho3sigma(npoints, 12);
                CDenseMatrix v4rho3lapl(npoints, 8);
                CDenseMatrix v4rho3tau(npoints, 8);
                CDenseMatrix v4rho2sigma2(npoints, 18);
                CDenseMatrix v4rho2sigmalapl(npoints, 18);
                CDenseMatrix v4rho2sigmatau(npoints, 18);
                CDenseMatrix v4rho2lapl2(npoints, 9);
                CDenseMatrix v4rho2lapltau(npoints, 12);
                CDenseMatrix v4rho2tau2(npoints, 9);
                CDenseMatrix v4rhosigma3(npoints, 20);
                // v4rhosigma2lapl: inconsistent size in libxc (36 vs 24);
                CDenseMatrix v4rhosigma2lapl(npoints, 36);
                // v4rhosigma2tau: inconsistent size in libxc (36 vs 24);
                CDenseMatrix v4rhosigma2tau(npoints, 36);
                CDenseMatrix v4rhosigmalapl2(npoints, 18);
                CDenseMatrix v4rhosigmalapltau(npoints, 24);
                // v4rhosigmatau2: inconsistent size in libxc (36 vs 18);
                CDenseMatrix v4rhosigmatau2(npoints, 36);
                CDenseMatrix v4rholapl3(npoints, 8);
                CDenseMatrix v4rholapl2tau(npoints, 12);
                CDenseMatrix v4rholapltau2(npoints, 12);
                CDenseMatrix v4rhotau3(npoints, 8);
                CDenseMatrix v4sigma4(npoints, 15);
                CDenseMatrix v4sigma3lapl(npoints, 20);
                // v4sigma3tau: inconsistent size in libxc (30 vs 20);
                CDenseMatrix v4sigma3tau(npoints, 30);
                CDenseMatrix v4sigma2lapl2(npoints, 18);
                CDenseMatrix v4sigma2lapltau(npoints, 24);
                CDenseMatrix v4sigma2tau2(npoints, 18);
                CDenseMatrix v4sigmalapl3(npoints, 12);
                CDenseMatrix v4sigmalapl2tau(npoints, 18);
                CDenseMatrix v4sigmalapltau2(npoints, 18);
                CDenseMatrix v4sigmatau3(npoints, 12);
                CDenseMatrix v4lapl4(npoints, 5);
                CDenseMatrix v4lapl3tau(npoints, 8);
                CDenseMatrix v4lapl2tau2(npoints, 9);
                CDenseMatrix v4lapltau3(npoints, 8);
                CDenseMatrix v4tau4(npoints, 5);
                auto         func = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                func.compute_lxc_for_mgga(npoints,
                                          rho.data(),
                                          sigma.data(),
                                          lapl.data(),
                                          tau.data(),
                                          v4rho4.values(),
                                          v4rho3sigma.values(),
                                          v4rho3lapl.values(),
                                          v4rho3tau.values(),
                                          v4rho2sigma2.values(),
                                          v4rho2sigmalapl.values(),
                                          v4rho2sigmatau.values(),
                                          v4rho2lapl2.values(),
                                          v4rho2lapltau.values(),
                                          v4rho2tau2.values(),
                                          v4rhosigma3.values(),
                                          v4rhosigma2lapl.values(),
                                          v4rhosigma2tau.values(),
                                          v4rhosigmalapl2.values(),
                                          v4rhosigmalapltau.values(),
                                          v4rhosigmatau2.values(),
                                          v4rholapl3.values(),
                                          v4rholapl2tau.values(),
                                          v4rholapltau2.values(),
                                          v4rhotau3.values(),
                                          v4sigma4.values(),
                                          v4sigma3lapl.values(),
                                          v4sigma3tau.values(),
                                          v4sigma2lapl2.values(),
                                          v4sigma2lapltau.values(),
                                          v4sigma2tau2.values(),
                                          v4sigmalapl3.values(),
                                          v4sigmalapl2tau.values(),
                                          v4sigmalapltau2.values(),
                                          v4sigmatau3.values(),
                                          v4lapl4.values(),
                                          v4lapl3tau.values(),
                                          v4lapl2tau2.values(),
                                          v4lapltau3.values(),
                                          v4tau4.values());
                py::dict ret;
                ret["v4rho4"] = vlx_general::pointer_to_numpy(v4rho4.values(), v4rho4.getNumberOfRows(), v4rho4.getNumberOfColumns());
                ret["v4rho3sigma"] =
                    vlx_general::pointer_to_numpy(v4rho3sigma.values(), v4rho3sigma.getNumberOfRows(), v4rho3sigma.getNumberOfColumns());
                ret["v4rho3lapl"] = vlx_general::pointer_to_numpy(v4rho3lapl.values(), v4rho3lapl.getNumberOfRows(), v4rho3lapl.getNumberOfColumns());
                ret["v4rho3tau"]  = vlx_general::pointer_to_numpy(v4rho3tau.values(), v4rho3tau.getNumberOfRows(), v4rho3tau.getNumberOfColumns());
                ret["v4rho2sigma2"] =
                    vlx_general::pointer_to_numpy(v4rho2sigma2.values(), v4rho2sigma2.getNumberOfRows(), v4rho2sigma2.getNumberOfColumns());
                ret["v4rho2sigmalapl"] =
                    vlx_general::pointer_to_numpy(v4rho2sigmalapl.values(), v4rho2sigmalapl.getNumberOfRows(), v4rho2sigmalapl.getNumberOfColumns());
                ret["v4rho2sigmatau"] =
                    vlx_general::pointer_to_numpy(v4rho2sigmatau.values(), v4rho2sigmatau.getNumberOfRows(), v4rho2sigmatau.getNumberOfColumns());
                ret["v4rho2lapl2"] =
                    vlx_general::pointer_to_numpy(v4rho2lapl2.values(), v4rho2lapl2.getNumberOfRows(), v4rho2lapl2.getNumberOfColumns());
                ret["v4rho2lapltau"] =
                    vlx_general::pointer_to_numpy(v4rho2lapltau.values(), v4rho2lapltau.getNumberOfRows(), v4rho2lapltau.getNumberOfColumns());
                ret["v4rho2tau2"] = vlx_general::pointer_to_numpy(v4rho2tau2.values(), v4rho2tau2.getNumberOfRows(), v4rho2tau2.getNumberOfColumns());
                ret["v4rhosigma3"] =
                    vlx_general::pointer_to_numpy(v4rhosigma3.values(), v4rhosigma3.getNumberOfRows(), v4rhosigma3.getNumberOfColumns());
                ret["v4rhosigma2lapl"] =
                    vlx_general::pointer_to_numpy(v4rhosigma2lapl.values(), v4rhosigma2lapl.getNumberOfRows(), v4rhosigma2lapl.getNumberOfColumns());
                ret["v4rhosigma2tau"] =
                    vlx_general::pointer_to_numpy(v4rhosigma2tau.values(), v4rhosigma2tau.getNumberOfRows(), v4rhosigma2tau.getNumberOfColumns());
                ret["v4rhosigmalapl2"] =
                    vlx_general::pointer_to_numpy(v4rhosigmalapl2.values(), v4rhosigmalapl2.getNumberOfRows(), v4rhosigmalapl2.getNumberOfColumns());
                ret["v4rhosigmalapltau"] = vlx_general::pointer_to_numpy(
                    v4rhosigmalapltau.values(), v4rhosigmalapltau.getNumberOfRows(), v4rhosigmalapltau.getNumberOfColumns());
                ret["v4rhosigmatau2"] =
                    vlx_general::pointer_to_numpy(v4rhosigmatau2.values(), v4rhosigmatau2.getNumberOfRows(), v4rhosigmatau2.getNumberOfColumns());
                ret["v4rholapl3"] = vlx_general::pointer_to_numpy(v4rholapl3.values(), v4rholapl3.getNumberOfRows(), v4rholapl3.getNumberOfColumns());
                ret["v4rholapl2tau"] =
                    vlx_general::pointer_to_numpy(v4rholapl2tau.values(), v4rholapl2tau.getNumberOfRows(), v4rholapl2tau.getNumberOfColumns());
                ret["v4rholapltau2"] =
                    vlx_general::pointer_to_numpy(v4rholapltau2.values(), v4rholapltau2.getNumberOfRows(), v4rholapltau2.getNumberOfColumns());
                ret["v4rhotau3"] = vlx_general::pointer_to_numpy(v4rhotau3.values(), v4rhotau3.getNumberOfRows(), v4rhotau3.getNumberOfColumns());
                ret["v4sigma4"]  = vlx_general::pointer_to_numpy(v4sigma4.values(), v4sigma4.getNumberOfRows(), v4sigma4.getNumberOfColumns());
                ret["v4sigma3lapl"] =
                    vlx_general::pointer_to_numpy(v4sigma3lapl.values(), v4sigma3lapl.getNumberOfRows(), v4sigma3lapl.getNumberOfColumns());
                ret["v4sigma3tau"] =
                    vlx_general::pointer_to_numpy(v4sigma3tau.values(), v4sigma3tau.getNumberOfRows(), v4sigma3tau.getNumberOfColumns());
                ret["v4sigma2lapl2"] =
                    vlx_general::pointer_to_numpy(v4sigma2lapl2.values(), v4sigma2lapl2.getNumberOfRows(), v4sigma2lapl2.getNumberOfColumns());
                ret["v4sigma2lapltau"] =
                    vlx_general::pointer_to_numpy(v4sigma2lapltau.values(), v4sigma2lapltau.getNumberOfRows(), v4sigma2lapltau.getNumberOfColumns());
                ret["v4sigma2tau2"] =
                    vlx_general::pointer_to_numpy(v4sigma2tau2.values(), v4sigma2tau2.getNumberOfRows(), v4sigma2tau2.getNumberOfColumns());
                ret["v4sigmalapl3"] =
                    vlx_general::pointer_to_numpy(v4sigmalapl3.values(), v4sigmalapl3.getNumberOfRows(), v4sigmalapl3.getNumberOfColumns());
                ret["v4sigmalapl2tau"] =
                    vlx_general::pointer_to_numpy(v4sigmalapl2tau.values(), v4sigmalapl2tau.getNumberOfRows(), v4sigmalapl2tau.getNumberOfColumns());
                ret["v4sigmalapltau2"] =
                    vlx_general::pointer_to_numpy(v4sigmalapltau2.values(), v4sigmalapltau2.getNumberOfRows(), v4sigmalapltau2.getNumberOfColumns());
                ret["v4sigmatau3"] =
                    vlx_general::pointer_to_numpy(v4sigmatau3.values(), v4sigmatau3.getNumberOfRows(), v4sigmatau3.getNumberOfColumns());
                ret["v4lapl4"]    = vlx_general::pointer_to_numpy(v4lapl4.values(), v4lapl4.getNumberOfRows(), v4lapl4.getNumberOfColumns());
                ret["v4lapl3tau"] = vlx_general::pointer_to_numpy(v4lapl3tau.values(), v4lapl3tau.getNumberOfRows(), v4lapl3tau.getNumberOfColumns());
                ret["v4lapl2tau2"] =
                    vlx_general::pointer_to_numpy(v4lapl2tau2.values(), v4lapl2tau2.getNumberOfRows(), v4lapl2tau2.getNumberOfColumns());
                ret["v4lapltau3"] = vlx_general::pointer_to_numpy(v4lapltau3.values(), v4lapltau3.getNumberOfRows(), v4lapltau3.getNumberOfColumns());
                ret["v4tau4"]     = vlx_general::pointer_to_numpy(v4tau4.values(), v4tau4.getNumberOfRows(), v4tau4.getNumberOfColumns());
                return ret;
            },
            "Computes Lxc for meta-GGA.",
            "xcFuncLabel"_a,
            "rho"_a,
            "sigma"_a,
            "lapl"_a,
            "tau"_a);

    // CXCMolecularGradient class

    PyClass<CXCMolecularGradient>(m, "XCMolecularGradient")
        .def(py::init(&vlx_general::create<CXCMolecularGradient>), "comm"_a = py::none())
        .def(
            "integrate_vxc_gradient",
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
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
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& rwDensityMatrix,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
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
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& rwDensityMatrixOne,
               const CAODensityMatrix& rwDensityMatrixTwo,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
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
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& rwDensityMatrixOne,
               const CAODensityMatrix& rwDensityMatrixTwo,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
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

    // CXCMolecularHessian class

    PyClass<CXCMolecularHessian>(m, "XCMolecularHessian")
        .def(py::init(&vlx_general::create<CXCMolecularHessian>), "comm"_a = py::none())
        .def(
            "integrate_exc_hessian",
            [](CXCMolecularHessian&    self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
                auto molhess = self.integrateExcHessian(molecule, basis, gsDensityMatrix, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molhess.values(), molhess.getNumberOfRows(), molhess.getNumberOfColumns());
            },
            "Integrates XC contribution to molecular Hessian.",
            "molecule"_a,
            "basis"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_vxc_fock_gradient",
            [](CXCMolecularHessian&    self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CAODensityMatrix& gsDensityMatrix,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel,
               const int32_t           atomIdx) -> py::array_t<double> {
                auto     vxcgrad = self.integrateVxcFockGradient(molecule, basis, gsDensityMatrix, molecularGrid, xcFuncLabel, atomIdx);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(vxcgrad[0].values(), vxcgrad[0].getNumberOfRows(), vxcgrad[0].getNumberOfColumns()));
                ret.append(vlx_general::pointer_to_numpy(vxcgrad[1].values(), vxcgrad[1].getNumberOfRows(), vxcgrad[1].getNumberOfColumns()));
                ret.append(vlx_general::pointer_to_numpy(vxcgrad[2].values(), vxcgrad[2].getNumberOfRows(), vxcgrad[2].getNumberOfColumns()));
                return ret;
            },
            "Integrates XC contribution to gradient of Vxc matrix element.",
            "molecule"_a,
            "basis"_a,
            "gsDensityMatrix"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a,
            "atomIdx"_a);

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

    m.def("available_functionals", &vxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
