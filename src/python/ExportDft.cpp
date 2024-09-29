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

// constructors for CGridDriver, CXCIntegrator

static auto
CGridDriver_create(py::object py_comm) -> std::shared_ptr<CGridDriver>
{
    if (py_comm.is_none()) return std::make_shared<CGridDriver>(MPI_COMM_WORLD);

    return std::make_shared<CGridDriver>(*vlx_general::get_mpi_comm(py_comm));
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
        .def(py::init<const CMolecularGrid&>())
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
        .def(py::self == py::self);

    // CGridDriver class

    PyClass<CGridDriver>(m, "GridDriver")
        .def(py::init(&CGridDriver_create), "comm"_a = py::none())
        .def("generate", &CGridDriver::generate, "Generates molecular grid for molecule.", "molecule"_a)
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
                auto fock_pointers    = arrays_to_mutable_pointers(aoFockArrays);
                auto rw_dens_pointers = arrays_to_const_pointers(rwDensityArrays);
                auto gs_dens_pointers = arrays_to_const_pointers(gsDensityArrays);
                self.integrateFxcFock(fock_pointers, molecule, basis, rw_dens_pointers, gs_dens_pointers, molecularGrid, fvxc);
            },
            "Integrates 2nd-order exchange-correlation contribution.")
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

    // exposing functions

    m.def("available_functionals", &vxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
