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

    for (size_t i = 0; i < arrays.size(); i++) {
        auto c_style = py::detail::check_flags(arrays[i].ptr(), py::array::c_style);
        errors::assertMsgCritical(c_style, errstyle);

        errors::assertMsgCritical(static_cast<int>(arrays[i].shape(0)) == nao, errshape);
        errors::assertMsgCritical(static_cast<int>(arrays[i].shape(1)) == nao, errshape);
    }
}

static auto
integrate_vxc_pdft(const CXCIntegrator&       self,
                   const CAODensityMatrix&    aoDensityMatrix,
                   const py::array_t<double>& active2DM,
                   const py::array_t<double>& activeMOs,
                   const CMolecule&           molecule,
                   const CMolecularBasis&     basis,
                   const CMolecularGrid&      molecularGrid,
                   const std::string&         xcFuncLabel) -> py::list
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

    CDenseMatrix tensor2DM(n_active * n_active, n_active * n_active);

    std::memcpy(tensor2DM.values(), active2DM.data(), active2DM.size() * sizeof(double));

    // activeMOs

    // Check dimensions

    errdim = "integrate_vxc_pdft, activeMOs: Expecting a 2D numpy array";

    errors::assertMsgCritical(activeMOs.ndim() == 2, errdim);

    // check that the numpy array is c-style contiguous

    errsrc = "integrate_vxc_pdft, activeMOs: Expecting a C-style contiguous numpy array";

    c_style = py::detail::check_flags(activeMOs.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    auto naos = activeMOs.shape(1);

    CDenseMatrix denseActiveMO(n_active, naos);

    std::memcpy(denseActiveMO.values(), activeMOs.data(), activeMOs.size() * sizeof(double));

    // Create output tensors

    CAOKohnShamMatrix matrixVxc(naos, naos, true);

    matrixVxc.zero();

    CDense4DTensor tensorWxc(naos, n_active, n_active, n_active);

    tensorWxc.zero();

    self.integrateVxcPDFT(matrixVxc, tensorWxc, molecule, basis, aoDensityMatrix, tensor2DM, denseActiveMO, molecularGrid, xcFuncLabel);

    py::list returnList;

    returnList.append(matrixVxc);

    returnList.append(vlx_general::pointer_to_numpy(tensorWxc.values(), {naos, n_active * n_active * n_active}));

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
        .def(
            "_generate_local_grid",
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
                auto numdensities = static_cast<int>(gsDensityArrays.size());
                std::string errsize("integrate_vxc_fock: Expecting a list of 1 or 2 numpy arrays");
                errors::assertMsgCritical((numdensities == 1) || (numdensities == 2), errsize);
                auto nao = basis.dimensions_of_basis();
                check_arrays("integrate_vxc_fock", gsDensityArrays, nao);
                auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
                return self.integrateVxcFock(molecule, basis, gs_dens_pointers, molecularGrid, xcFuncLabel);
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
                auto numdensities = static_cast<int>(gsDensityArrays.size());
                std::string errsize("integrate_vxc_fock: Expecting a list of 1 or 2 numpy arrays");
                errors::assertMsgCritical((numdensities == 1) || (numdensities == 2), errsize);
                auto nao = basis.dimensions_of_basis();
                check_arrays("integrate_vxc_fock", gsDensityArrays, nao);
                auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
                return self.integrateVxcFock(molecule, basis, gs_dens_pointers, molecularGrid, fvxc);
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
                auto num_focks = static_cast<int>(aoFockArrays.size());
                auto num_rw_dens = static_cast<int>(rwDensityArrays.size());
                auto num_gs_dens = static_cast<int>(gsDensityArrays.size());
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
                self.integrateFxcFock(fock_pointers, molecule, basis, rw_dens_pointers, gs_dens_pointers, molecularGrid, xcFuncLabel);
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
                auto num_focks = static_cast<int>(aoFockArrays.size());
                auto num_rw_dens = static_cast<int>(rwDensityArrays.size());
                auto num_gs_dens = static_cast<int>(gsDensityArrays.size());
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
            },
            "Integrates 2nd-order exchange-correlation contribution.")
        .def("integrate_vxc_pdft", &integrate_vxc_pdft)
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
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const py::array_t<double>& gsDensity,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
                auto gsDensityPointer = gsDensity.data();
                auto molgrad = self.integrateVxcGradient(molecule, basis, gsDensityPointer, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), {molgrad.getNumberOfRows(), molgrad.getNumberOfColumns()});
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "gsDensity"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a)
        .def(
            "integrate_vxc_gradient",
            [](CXCMolecularGradient&   self,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const py::array_t<double>& rwDensity,
               const py::array_t<double>& gsDensity,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFuncLabel) -> py::array_t<double> {
                auto rwDensityPointer = rwDensity.data();
                auto gsDensityPointer = gsDensity.data();
                auto molgrad = self.integrateVxcGradient(molecule, basis, rwDensityPointer, gsDensityPointer, molecularGrid, xcFuncLabel);
                return vlx_general::pointer_to_numpy(molgrad.values(), {molgrad.getNumberOfRows(), molgrad.getNumberOfColumns()});
            },
            "Integrates 1st-order exchange-correlation contribution to molecular gradient.",
            "molecule"_a,
            "basis"_a,
            "rwDensity"_a,
            "gsDensity"_a,
            "molecularGrid"_a,
            "xcFuncLabel"_a);

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

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
