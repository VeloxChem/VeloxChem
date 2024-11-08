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

#include "ExportDft.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <vector>

#include "AOKohnShamMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "FunctionalParser.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCComponent.hpp"
#include "XCFunctional.hpp"
#include "XCIntegrator.hpp"
#include "XCMolecularGradient.hpp"
#include "XCMolecularHessian.hpp"
#include "XCPairDensityFunctional.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

// helper functions for converting arrays to pointers

static auto
arrays_to_mutable_pointers(std::vector<py::array_t<double>>& arrays) -> std::vector<double*>
{
    std::vector<double*> pointers;
    for (size_t i = 0; i < arrays.size(); i++)
    {
        std::string errstyle("arrays_to_mutable_pointers: Expecting contiguous numpy arrays");
        auto        c_style = py::detail::check_flags(arrays[i].ptr(), py::array::c_style);
        errors::assertMsgCritical(c_style, errstyle);
        pointers.push_back(arrays[i].mutable_data());
    }
    return pointers;
}

static auto
arrays_to_const_pointers(const std::vector<py::array_t<double>>& arrays) -> std::vector<const double*>
{
    std::vector<const double*> pointers;
    for (size_t i = 0; i < arrays.size(); i++)
    {
        std::string errstyle("arrays_to_const_pointers: Expecting contiguous numpy arrays");
        auto        c_style = py::detail::check_flags(arrays[i].ptr(), py::array::c_style);
        errors::assertMsgCritical(c_style, errstyle);
        pointers.push_back(arrays[i].data());
    }
    return pointers;
}

static auto
check_arrays(const std::string& func_name, const std::vector<py::array_t<double>>& arrays, const int nao) -> void
{
    std::string errstyle(func_name + std::string(": Expecting contiguous numpy arrays"));
    std::string errshape(func_name + std::string(": Invalide shape of numpy array"));

    for (size_t i = 0; i < arrays.size(); i++)
    {
        auto c_style = py::detail::check_flags(arrays[i].ptr(), py::array::c_style);
        errors::assertMsgCritical(c_style, errstyle);

        errors::assertMsgCritical(static_cast<int>(arrays[i].shape(0)) == nao, errshape);
        errors::assertMsgCritical(static_cast<int>(arrays[i].shape(1)) == nao, errshape);
    }
}

static auto
CXCIntegrator_integrate_vxc_fock(const CXCIntegrator&                    self,
                                 const CMolecule&                        molecule,
                                 const CMolecularBasis&                  basis,
                                 const std::vector<py::array_t<double>>& gsDensityArrays,
                                 const CMolecularGrid&                   molecularGrid,
                                 const CXCFunctional&                    fvxc) -> CAOKohnShamMatrix
{
    auto        numdensities = static_cast<int>(gsDensityArrays.size());
    std::string errsize("integrate_vxc_fock: Expecting a list of 1 or 2 numpy arrays");
    errors::assertMsgCritical((numdensities == 1) || (numdensities == 2), errsize);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_vxc_fock", gsDensityArrays, nao);

    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    return self.integrateVxcFock(molecule, basis, gs_dens_pointers, molecularGrid, fvxc);
}

static auto
CXCIntegrator_integrate_fxc_fock(const CXCIntegrator&                    self,
                                 std::vector<py::array_t<double>>&       aoFockArrays,
                                 const CMolecule&                        molecule,
                                 const CMolecularBasis&                  basis,
                                 const std::vector<py::array_t<double>>& rwDensityArrays,
                                 const std::vector<py::array_t<double>>& gsDensityArrays,
                                 const CMolecularGrid&                   molecularGrid,
                                 const CXCFunctional&                    fvxc) -> void
{
    auto        num_focks   = static_cast<int>(aoFockArrays.size());
    auto        num_rw_dens = static_cast<int>(rwDensityArrays.size());
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_fxc_fock: Inconsistent number of numpy arrays");
    errors::assertMsgCritical(num_rw_dens == num_focks, errnum);
    errors::assertMsgCritical(num_gs_dens == 1, errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_fxc_fock", aoFockArrays, nao);
    check_arrays("integrate_fxc_fock", rwDensityArrays, nao);
    check_arrays("integrate_fxc_fock", gsDensityArrays, nao);

    auto fock_pointers    = arrays_to_mutable_pointers(aoFockArrays);
    auto rw_dens_pointers = arrays_to_const_pointers(rwDensityArrays);
    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    self.integrateFxcFock(fock_pointers, molecule, basis, rw_dens_pointers, gs_dens_pointers, molecularGrid, fvxc);
}

static auto
CXCIntegrator_integrate_kxc_fock(const CXCIntegrator&                    self,
                                 std::vector<py::array_t<double>>&       aoFockArrays,
                                 const CMolecule&                        molecule,
                                 const CMolecularBasis&                  basis,
                                 const std::vector<py::array_t<double>>& rwDensityArrays,
                                 const std::vector<py::array_t<double>>& rw2DensityArrays,
                                 const std::vector<py::array_t<double>>& gsDensityArrays,
                                 const CMolecularGrid&                   molecularGrid,
                                 const CXCFunctional&                    fvxc,
                                 const std::string&                      quadMode) -> void
{
    auto        num_focks    = static_cast<int>(aoFockArrays.size());
    auto        num_rw2_dens = static_cast<int>(rw2DensityArrays.size());
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_kxc_fock: Inconsistent number of numpy arrays");
    errors::assertMsgCritical(num_rw2_dens == num_focks, errnum);
    errors::assertMsgCritical(num_gs_dens == 1, errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_kxc_fock", aoFockArrays, nao);
    check_arrays("integrate_kxc_fock", rwDensityArrays, nao);
    check_arrays("integrate_kxc_fock", rw2DensityArrays, nao);
    check_arrays("integrate_kxc_fock", gsDensityArrays, nao);

    auto fock_pointers = arrays_to_mutable_pointers(aoFockArrays);
    auto rw_dens_pointers = arrays_to_const_pointers(rwDensityArrays);
    auto rw2_dens_pointers = arrays_to_const_pointers(rw2DensityArrays);
    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    self.integrateKxcFock(fock_pointers, molecule, basis, rw_dens_pointers, rw2_dens_pointers, gs_dens_pointers, molecularGrid, fvxc, quadMode);
}

static auto
CXCIntegrator_integrate_kxclxc_fock(const CXCIntegrator&              self,
                                    std::vector<py::array_t<double>>& aoFockArrays,
                                    const CMolecule&                  molecule,
                                    const CMolecularBasis&            basis,
                                    const std::vector<py::array_t<double>>& rwDensityArrays,
                                    const std::vector<py::array_t<double>>& rw2DensityArrays,
                                    const std::vector<py::array_t<double>>& rw3DensityArrays,
                                    const std::vector<py::array_t<double>>& gsDensityArrays,
                                    const CMolecularGrid&             molecularGrid,
                                    const CXCFunctional&              fvxc,
                                    const std::string&                cubeMode) -> void
{
    auto        num_focks    = static_cast<int>(aoFockArrays.size());
    auto        num_rw2_dens = static_cast<int>(rw2DensityArrays.size());
    auto        num_rw3_dens = static_cast<int>(rw3DensityArrays.size());
    auto        num_gs_dens  = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_kxclxc_fock: Inconsistent number of numpy arrays");
    errors::assertMsgCritical(num_rw2_dens + num_rw3_dens == num_focks, errnum);
    errors::assertMsgCritical(num_gs_dens == 1, errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_kxclxc_fock", aoFockArrays, nao);
    check_arrays("integrate_kxclxc_fock", rwDensityArrays, nao);
    check_arrays("integrate_kxclxc_fock", rw2DensityArrays, nao);
    check_arrays("integrate_kxclxc_fock", rw3DensityArrays, nao);
    check_arrays("integrate_kxclxc_fock", gsDensityArrays, nao);

    auto fock_pointers = arrays_to_mutable_pointers(aoFockArrays);
    auto rw_dens_pointers = arrays_to_const_pointers(rwDensityArrays);
    auto rw2_dens_pointers = arrays_to_const_pointers(rw2DensityArrays);
    auto rw3_dens_pointers = arrays_to_const_pointers(rw3DensityArrays);
    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    self.integrateKxcLxcFock(fock_pointers, molecule, basis, rw_dens_pointers, rw2_dens_pointers, rw3_dens_pointers, gs_dens_pointers, molecularGrid, fvxc, cubeMode);
}

static auto
CXCIntegrator_integrate_vxc_pdft(const CXCIntegrator&       self,
                   const py::array_t<double>& densityMatrix,
                   const py::array_t<double>& active2DM,
                   const py::array_t<double>& activeMOs,
                   const CMolecule&           molecule,
                   const CMolecularBasis&     basis,
                   const CMolecularGrid&      molecularGrid,
                   const std::string&         xcFuncLabel,
                   const py::dict             components,
                   const double               rs_omega) -> py::list
{
    // active2DM

    // sanity checks
    std::string errsrc("integrate_vxc_pdft, active2DM: Expecting a C-style contiguous numpy array");
    auto        c_style = py::detail::check_flags(active2DM.ptr(), py::array::c_style);
    errors::assertMsgCritical(c_style, errsrc);

    std::string errsizes("integrate_vxc_pdft, active2DM: Expecting 4 identical dimensions");
    bool        same_size = ((active2DM.ndim() == 4) && (active2DM.shape(0) == active2DM.shape(1)) && (active2DM.shape(0) == active2DM.shape(2)) &&
                      (active2DM.shape(0) == active2DM.shape(3)));
    errors::assertMsgCritical(same_size, errsizes);

    auto n_active = static_cast<int>(active2DM.shape(0));

    // Form 4D tensor
    CDenseMatrix tensor2DM(n_active * n_active, n_active * n_active);
    std::memcpy(tensor2DM.values(), active2DM.data(), active2DM.size() * sizeof(double));

    // activeMOs

    // sanity checks
    errsrc  = "integrate_vxc_pdft, activeMOs: Expecting a C-style contiguous numpy array";
    c_style = py::detail::check_flags(activeMOs.ptr(), py::array::c_style);
    errors::assertMsgCritical(c_style, errsrc);

    errsizes = "integrate_vxc_pdft, activeMOs: Expecting a 2D numpy array";
    errors::assertMsgCritical(activeMOs.ndim() == 2, errsizes);

    auto naos = static_cast<int>(activeMOs.shape(1));

    CDenseMatrix denseActiveMO(n_active, naos);
    std::memcpy(denseActiveMO.values(), activeMOs.data(), activeMOs.size() * sizeof(double));

    // densityMatrix

    // sanity checks
    errsrc  = "integrate_vxc_pdft, densityMatrix: Expecting a C-style contiguous numpy array";
    c_style = py::detail::check_flags(densityMatrix.ptr(), py::array::c_style);
    errors::assertMsgCritical(c_style, errsrc);

    errsizes  = "integrate_vxc_pdft, densityMatrix: Expecting 2 identical dimensions";
    same_size = ((densityMatrix.ndim() == 2) && (densityMatrix.shape(0) == naos) && (densityMatrix.shape(1) == naos));
    errors::assertMsgCritical(same_size, errsizes);

    // Functional

    std::vector<std::string> labels;
    std::vector<double> coeffs;
    for (auto item : components)
    {
        labels.push_back(item.first.cast<std::string>());
        coeffs.push_back(item.second.cast<double>());
    };
    CXCPairDensityFunctional xcfun = CXCPairDensityFunctional(xcFuncLabel, labels, coeffs);

    // Create output tensors

    CAOKohnShamMatrix matrixVxc(naos, naos, true);
    matrixVxc.zero();

    CDense4DTensor tensorWxc(naos, n_active, n_active, n_active);
    tensorWxc.zero();

    self.integrateVxcPDFT(matrixVxc, tensorWxc, molecule, basis, densityMatrix.data(), tensor2DM, denseActiveMO, molecularGrid, xcfun, rs_omega);

    py::list returnList;
    returnList.append(matrixVxc);
    returnList.append(vlx_general::pointer_to_numpy(tensorWxc.values(), {naos, n_active, n_active, n_active}));
    return returnList;
}

auto
CXCIntegrator_integrate_vxc_gradient(CXCMolecularGradient&      self,
                                     const CMolecule&           molecule,
                                     const CMolecularBasis&     basis,
                                     const std::vector<py::array_t<double>>& rwDensityArrays,
                                     const std::vector<py::array_t<double>>& gsDensityArrays,
                                     const CMolecularGrid&      molecularGrid,
                                     const std::string&         xcFuncLabel) -> py::array_t<double>
{
    auto        num_rw_dens = static_cast<int>(rwDensityArrays.size());
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_vxc_gradient: Inconsistent number of numpy arrays");
    errors::assertMsgCritical(num_rw_dens == num_gs_dens, errnum);
    errors::assertMsgCritical((num_gs_dens == 1) || (num_gs_dens == 2), errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_vxc_gradient", rwDensityArrays, nao);
    check_arrays("integrate_vxc_gradient", gsDensityArrays, nao);

    auto rw_dens_pointers = arrays_to_const_pointers(rwDensityArrays);
    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    auto molgrad          = self.integrateVxcGradient(molecule, basis, rw_dens_pointers, gs_dens_pointers, molecularGrid, xcFuncLabel);
    return vlx_general::pointer_to_numpy(molgrad.values(), {molgrad.getNumberOfRows(), molgrad.getNumberOfColumns()});
}

auto
CXCIntegrator_integrate_fxc_gradient(CXCMolecularGradient&      self,
                                     const CMolecule&           molecule,
                                     const CMolecularBasis&     basis,
                                     const std::vector<py::array_t<double>>& rwDensityArraysOne,
                                     const std::vector<py::array_t<double>>& rwDensityArraysTwo,
                                     const std::vector<py::array_t<double>>& gsDensityArrays,
                                     const CMolecularGrid&      molecularGrid,
                                     const std::string&         xcFuncLabel) -> py::array_t<double>
{
    auto        num_rw_dens_one = static_cast<int>(rwDensityArraysOne.size());
    auto        num_rw_dens_two = static_cast<int>(rwDensityArraysTwo.size());
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_fxc_gradient: Inconsistent number of numpy arrays");
    errors::assertMsgCritical((num_rw_dens_one == num_gs_dens) && (num_rw_dens_two == num_gs_dens), errnum);
    errors::assertMsgCritical((num_gs_dens == 1) || (num_gs_dens == 2), errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_fxc_gradient", rwDensityArraysOne, nao);
    check_arrays("integrate_fxc_gradient", rwDensityArraysTwo, nao);
    check_arrays("integrate_fxc_gradient", gsDensityArrays, nao);

    auto rw_dens_pointers_one = arrays_to_const_pointers(rwDensityArraysOne);
    auto rw_dens_pointers_two = arrays_to_const_pointers(rwDensityArraysTwo);
    auto gs_dens_pointers     = arrays_to_const_pointers(gsDensityArrays);
    auto molgrad              = self.integrateFxcGradient(molecule, basis, rw_dens_pointers_one, rw_dens_pointers_two, gs_dens_pointers, molecularGrid, xcFuncLabel);
    return vlx_general::pointer_to_numpy(molgrad.values(), {molgrad.getNumberOfRows(), molgrad.getNumberOfColumns()});
}

auto
CXCIntegrator_integrate_kxc_gradient(CXCMolecularGradient&      self,
                                     const CMolecule&           molecule,
                                     const CMolecularBasis&     basis,
                                     const std::vector<py::array_t<double>>& rwDensityArraysOne,
                                     const std::vector<py::array_t<double>>& rwDensityArraysTwo,
                                     const std::vector<py::array_t<double>>& gsDensityArrays,
                                     const CMolecularGrid&      molecularGrid,
                                     const std::string&         xcFuncLabel) -> py::array_t<double>
{
    auto        num_rw_dens_one = static_cast<int>(rwDensityArraysOne.size());
    auto        num_rw_dens_two = static_cast<int>(rwDensityArraysTwo.size());
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_kxc_gradient: Inconsistent number of numpy arrays");
    errors::assertMsgCritical((num_rw_dens_one == num_gs_dens) && (num_rw_dens_two == num_gs_dens), errnum);
    errors::assertMsgCritical((num_gs_dens == 1) || (num_gs_dens == 2), errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_kxc_gradient", rwDensityArraysOne, nao);
    check_arrays("integrate_kxc_gradient", rwDensityArraysTwo, nao);
    check_arrays("integrate_kxc_gradient", gsDensityArrays, nao);

    auto rw_dens_pointers_one = arrays_to_const_pointers(rwDensityArraysOne);
    auto rw_dens_pointers_two = arrays_to_const_pointers(rwDensityArraysTwo);
    auto gs_dens_pointers     = arrays_to_const_pointers(gsDensityArrays);
    auto molgrad              = self.integrateKxcGradient(molecule, basis, rw_dens_pointers_one, rw_dens_pointers_two, gs_dens_pointers, molecularGrid, xcFuncLabel);
    return vlx_general::pointer_to_numpy(molgrad.values(), {molgrad.getNumberOfRows(), molgrad.getNumberOfColumns()});
}

auto
CXCMolecularHessian_integrate_exc_hessian(CXCMolecularHessian&       self,
                                          const CMolecule&           molecule,
                                          const CMolecularBasis&     basis,
                                          const std::vector<py::array_t<double>>& gsDensityArrays,
                                          const CMolecularGrid&      molecularGrid,
                                          const std::string&         xcFuncLabel) -> py::array_t<double>
{
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_exc_hessian: Inconsistent number of numpy arrays");
    errors::assertMsgCritical((num_gs_dens == 1) || (num_gs_dens == 2), errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_exc_hessian", gsDensityArrays, nao);

    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    auto molhess          = self.integrateExcHessian(molecule, basis, gs_dens_pointers, molecularGrid, xcFuncLabel);
    return vlx_general::pointer_to_numpy(molhess.values(), {molhess.getNumberOfRows(), molhess.getNumberOfColumns()});
}

auto
CXCMolecularHessian_integrate_vxc_fock_gradient(CXCMolecularHessian&       self,
                                                const CMolecule&           molecule,
                                                const CMolecularBasis&     basis,
                                                const std::vector<py::array_t<double>>& gsDensityArrays,
                                                const CMolecularGrid&      molecularGrid,
                                                const std::string&         xcFuncLabel,
                                                const int                  atomIdx) -> py::array_t<double>
{
    auto        num_gs_dens = static_cast<int>(gsDensityArrays.size());
    std::string errnum("integrate_vxc_fock_gradient: Inconsistent number of numpy arrays");
    errors::assertMsgCritical((num_gs_dens == 1) || (num_gs_dens == 2), errnum);

    auto nao = basis.dimensions_of_basis();
    check_arrays("integrate_vxc_fock_gradient", gsDensityArrays, nao);

    auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
    auto vxcgrad          = self.integrateVxcFockGradient(molecule, basis, gs_dens_pointers, molecularGrid, xcFuncLabel, atomIdx);

    py::list ret;
    ret.append(vlx_general::pointer_to_numpy(vxcgrad[0].values(), {vxcgrad[0].getNumberOfRows(), vxcgrad[0].getNumberOfColumns()}));
    ret.append(vlx_general::pointer_to_numpy(vxcgrad[1].values(), {vxcgrad[1].getNumberOfRows(), vxcgrad[1].getNumberOfColumns()}));
    ret.append(vlx_general::pointer_to_numpy(vxcgrad[2].values(), {vxcgrad[2].getNumberOfRows(), vxcgrad[2].getNumberOfColumns()}));
    return ret;
}

// Exports classes/functions in src/dft to python

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
        .def(py::init<int, int, bool>(), "nrows"_a, "ncols"_a, "is_rest"_a)
        .def(
            "alpha_to_numpy",
            [](const CAOKohnShamMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaValues(), {self.getNumberOfRows(), self.getNumberOfColumns()});
            },
            "Converts alpha AOKohnShamMatrix to numpy array.")
        .def(
            "beta_to_numpy",
            [](const CAOKohnShamMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaValues(), {self.getNumberOfRows(), self.getNumberOfColumns()});
            },
            "Converts beta AOKohnShamMatrix to numpy array.")
        .def("get_electrons", &CAOKohnShamMatrix::getNumberOfElectrons, "Gets number of electrons obtained by integrating Kohn-Sham matrix.")
        .def("get_energy", &CAOKohnShamMatrix::getExchangeCorrelationEnergy, "Gets exchange-correlation energy associated with Kohn-Sham matrix.")
        .def(py::self == py::self);

    // CMolecularGrid class

    PyClass<CMolecularGrid>(m, "MolecularGrid")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init<const CMolecularGrid&>())
        .def("partition_grid_points", &CMolecularGrid::partitionGridPoints)
        .def("distribute_counts_and_displacements", &CMolecularGrid::distributeCountsAndDisplacements)
        .def("number_of_points", &CMolecularGrid::getNumberOfGridPoints)
        .def(
            "x_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesX(), std::vector<int>{self.getNumberOfGridPoints()});
            },
            "Gets X coordinates of grid as numpy array.")
        .def(
            "y_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesY(), std::vector<int>{self.getNumberOfGridPoints()});
            },
            "Gets Y coordinates of grid as numpy array.")
        .def(
            "z_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), std::vector<int>{self.getNumberOfGridPoints()});
            },
            "Gets Z coordinates of grid as numpy array.")
        .def(
            "w_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getWeights(), std::vector<int>{self.getNumberOfGridPoints()});
            },
            "Gets weights of grid as numpy array.")
        .def(
            "grid_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                auto points = self.getGridPoints();
                return vlx_general::pointer_to_numpy(points.values(), {4, self.getNumberOfGridPoints()});
            },
            "Gets grid points as numpy array of shape (4,N).")
        .def(py::self == py::self);

    // CGridDriver class
    // Note: GridDriver is prefixed by an underscore and will be used in griddriver.py

    PyClass<CGridDriver>(m, "_GridDriver")
        .def(py::init<>())
        .def("_generate_local_grid",
             &CGridDriver::generate_local_grid,
             "Generates MPI-local molecular grid for molecule.",
             "molecule"_a,
             "rank"_a,
             "nnodes"_a)
        .def("set_level", &CGridDriver::setLevel, "Sets accuracy level for grid generation.", "grid_level"_a);

    // CXCIntegrator class

    PyClass<CXCIntegrator>(m, "XCIntegrator")
        .def(py::init<>())
        .def(
            "integrate_vxc_fock",
            [](const CXCIntegrator&                    self,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const std::string&                      xcFuncLabel) -> CAOKohnShamMatrix {
                auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                return CXCIntegrator_integrate_vxc_fock(self, molecule, basis, gsDensityArrays, molecularGrid, fvxc);
            },
            "Integrates 1st-order exchange-correlation contribution.")
        .def(
            "integrate_vxc_fock",
            [](const CXCIntegrator&                    self,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const CXCFunctional&                    fvxc) -> CAOKohnShamMatrix {
                return CXCIntegrator_integrate_vxc_fock(self, molecule, basis, gsDensityArrays, molecularGrid, fvxc);
            },
            "Integrates 1st-order exchange-correlation contribution.")
        .def(
            "integrate_fxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const std::string&                      xcFuncLabel) -> void {
                auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                CXCIntegrator_integrate_fxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, gsDensityArrays, molecularGrid, fvxc);
            },
            "Integrates 2nd-order exchange-correlation contribution.")
        .def(
            "integrate_fxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const CXCFunctional&                    fvxc) -> void {
                CXCIntegrator_integrate_fxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, gsDensityArrays, molecularGrid, fvxc);
            },
            "Integrates 2nd-order exchange-correlation contribution.")
        .def(
            "integrate_kxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& rw2DensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const std::string&                      xcFuncLabel,
               const std::string&                      quadMode) -> void {
                auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                CXCIntegrator_integrate_kxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, rw2DensityArrays, gsDensityArrays, molecularGrid, fvxc, quadMode);
            },
            "Integrates 3rd-order exchange-correlation contribution.")
        .def(
            "integrate_kxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& rw2DensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const CXCFunctional&                    fvxc,
               const std::string&                      quadMode) -> void {
                CXCIntegrator_integrate_kxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, rw2DensityArrays, gsDensityArrays, molecularGrid, fvxc, quadMode);
            },
            "Integrates 3rd-order exchange-correlation contribution.")
        .def(
            "integrate_kxclxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& rw2DensityArrays,
               const std::vector<py::array_t<double>>& rw3DensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const std::string&                      xcFuncLabel,
               const std::string&                      cubeMode) -> void {
                auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
                CXCIntegrator_integrate_kxclxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, rw2DensityArrays, rw3DensityArrays, gsDensityArrays, molecularGrid, fvxc, cubeMode);
            },
            "Integrates 4th-order exchange-correlation contribution.")
        .def(
            "integrate_kxclxc_fock",
            [](const CXCIntegrator&                    self,
               std::vector<py::array_t<double>>&       aoFockArrays,
               const CMolecule&                        molecule,
               const CMolecularBasis&                  basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& rw2DensityArrays,
               const std::vector<py::array_t<double>>& rw3DensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&                   molecularGrid,
               const CXCFunctional&                    fvxc,
               const std::string&                      cubeMode) -> void {
                CXCIntegrator_integrate_kxclxc_fock(self, aoFockArrays, molecule, basis, rwDensityArrays, rw2DensityArrays, rw3DensityArrays, gsDensityArrays, molecularGrid, fvxc, cubeMode);
            },
            "Integrates 4th-order exchange-correlation contribution.")
        .def("integrate_vxc_pdft", &CXCIntegrator_integrate_vxc_pdft)
        .def(
            "compute_gto_values",
            [](CXCIntegrator& self, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid)
                -> py::array_t<double> {
                auto gtovalues = self.computeGtoValuesOnGridPoints(molecule, basis, molecularGrid);
                return vlx_general::pointer_to_numpy(gtovalues.values(), {gtovalues.getNumberOfRows(), gtovalues.getNumberOfColumns()});
            },
            "Computes GTO values on grid points.",
            "molecule"_a,
            "basis"_a,
            "molecular_grid"_a);

    // CXCMolecularGradient class

    PyClass<CXCMolecularGradient>(m, "XCMolecularGradient")
        .def(py::init<>())
        .def(
            "integrate_vxc_gradient",
            [](CXCMolecularGradient&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel) -> py::array_t<double> {
                return CXCIntegrator_integrate_vxc_gradient(self, molecule, basis, gsDensityArrays, gsDensityArrays, molecularGrid, xcFuncLabel);
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_vxc_gradient",
            [](CXCMolecularGradient&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& rwDensityArrays,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel) -> py::array_t<double> {
                return CXCIntegrator_integrate_vxc_gradient(self, molecule, basis, rwDensityArrays, gsDensityArrays, molecularGrid, xcFuncLabel);
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityArrays"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_fxc_gradient",
            [](CXCMolecularGradient&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& rwDensityArraysOne,
               const std::vector<py::array_t<double>>& rwDensityArraysTwo,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel) -> py::array_t<double> {
                return CXCIntegrator_integrate_fxc_gradient(self, molecule, basis, rwDensityArraysOne, rwDensityArraysTwo, gsDensityArrays, molecularGrid, xcFuncLabel);
            },
            "Integrates 2nd-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityArraysOne"_a,
            "rwDensityArraysTwo"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_kxc_gradient",
            [](CXCMolecularGradient&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& rwDensityArraysOne,
               const std::vector<py::array_t<double>>& rwDensityArraysTwo,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel) -> py::array_t<double> {
                return CXCIntegrator_integrate_kxc_gradient(self, molecule, basis, rwDensityArraysOne, rwDensityArraysTwo, gsDensityArrays, molecularGrid, xcFuncLabel);
            },
            "Integrates 3rd-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensityArraysOne"_a,
            "rwDensityArraysTwo"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a);

    // CXCMolecularHessian class

    PyClass<CXCMolecularHessian>(m, "XCMolecularHessian")
        .def(py::init<>())
        .def(
            "integrate_exc_hessian",
            [](CXCMolecularHessian&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel) -> py::array_t<double> {
                return CXCMolecularHessian_integrate_exc_hessian(self, molecule, basis, gsDensityArrays, molecularGrid, xcFuncLabel);
            },
            "Integrates exchange-correlation contribution to molecular Hessian.",
            "molecule"_a,
            "basis"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_vxc_fock_gradient",
            [](CXCMolecularHessian&      self,
               const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<py::array_t<double>>& gsDensityArrays,
               const CMolecularGrid&      molecularGrid,
               const std::string&         xcFuncLabel,
               const int                  atomIdx) -> py::array_t<double> {
                return CXCMolecularHessian_integrate_vxc_fock_gradient(self, molecule, basis, gsDensityArrays, molecularGrid, xcFuncLabel, atomIdx);
            },
            "Integrates exchange-correlation contribution to Vxc gradient.",
            "molecule"_a,
            "basis"_a,
            "gsDensityArrays"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a,
            "atomIdx"_a);

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
        .def("get_libxc_version", &CXCFunctional::getLibxcVersion, "Gets Libxc version.")
        .def("get_libxc_reference", &CXCFunctional::getLibxcReference, "Gets Libxc reference.")
        .def("get_functional_reference", &CXCFunctional::getFunctionalReference, "Gets functional reference.")
        .def("is_range_separated", &CXCFunctional::isRangeSeparated, "Determines whether the XC function is range-separated.")
        .def("is_hybrid", &CXCFunctional::isHybrid, "Determines whether the XC functional is hybrid.")
        .def("is_undefined", &CXCFunctional::isUndefined, "Determines whether the XC function is undefined.")
        .def("get_func_type", &CXCFunctional::getFunctionalType, "Gets type of XC functional.")
        .def("get_func_label", &CXCFunctional::getFunctionalLabel, "Gets name of XC functional.")
        .def("get_frac_exact_exchange", &CXCFunctional::getFractionOfExactExchange, "Gets fraction of exact Hartree-Fock exchange in XC functional.")
        .def("get_rs_alpha", &CXCFunctional::getRangeSeparationParameterAlpha, "Gets range-separation parameter alpha.")
        .def("get_rs_beta", &CXCFunctional::getRangeSeparationParameterBeta, "Gets range-separation parameter beta.")
        .def("get_rs_omega", &CXCFunctional::getRangeSeparationParameterOmega, "Gets range-separation parameter omega.")
        .def("get_dimension_of_derivatives", &CXCFunctional::getDimensionOfDerivatives, "Gets dimension of derivatives.")
        .def("set_rs_omega", &CXCFunctional::setRangeSeparatedParameterOmega, "Sets range-separation parameter omega.");

    // XCPairDensityFunctional class
    PyClass<CXCPairDensityFunctional>(m, "XCPairDensityFunctional")
        .def(py::init<const std::string&, const std::vector<std::string>&, const std::vector<double>&>(),
             "name_of_functional"_a,
             "labels"_a,
             "coeffs"_a)
        .def(py::init<const CXCPairDensityFunctional&>())
        .def(py::self == py::self)
        .def("get_func_label", &CXCPairDensityFunctional::getFunctionalLabel, "Gets functional name.")
        .def("get_func_type", &CXCPairDensityFunctional::getFunctionalType, "Gets functional type.");

    // exposing functions

    m.def("available_functionals", &vxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");
    m.def("available_pdft_functionals", &vxcfuncs::getAvailablePairDensityFunctionals, "Gets a list of available pdft exchange-correlation functionals components.");

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
