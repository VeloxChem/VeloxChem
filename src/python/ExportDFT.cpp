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

#include "ExportDFT.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <vector>

#include "AOKohnShamMatrix.hpp"
#include "ExportGeneral.hpp"
#include "FunctionalParser.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCComponent.hpp"
#include "XCIntegrator.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

// constructors for CGridDriver, CXCIntegrator

static auto
CGridDriver_create(py::object py_comm) -> std::shared_ptr<CGridDriver>
{
    if (py_comm.is_none()) return std::make_shared<CGridDriver>(MPI_COMM_WORLD);

    return std::make_shared<CGridDriver>(*vlx_general::get_mpi_comm(py_comm));
}

static auto
CXCIntegrator_create(py::object py_comm) -> std::shared_ptr<CXCIntegrator>
{
    if (py_comm.is_none()) return std::make_shared<CXCIntegrator>(MPI_COMM_WORLD);

    return std::make_shared<CXCIntegrator>(*vlx_general::get_mpi_comm(py_comm));
}

// Exports classes/functions in src/dft to python

void
export_dft(py::module& m)
{
    // CAOKohnShamMatrix class

    PyClass<CAOKohnShamMatrix>(m, "AOKohnShamMatrix")
        .def(py::init<>())
        .def(py::init<int64_t, int64_t, bool>(), "nrows"_a, "ncols"_a, "is_rest"_a)
        .def("get_alpha_matrix", &CAOKohnShamMatrix::getReferenceToAlphaKohnSham, "Gets constant reference to alpha-spin Kohn-Sham matrix.")
        .def("get_beta_matrix", &CAOKohnShamMatrix::getReferenceToBetaKohnSham, "Gets constant reference to beta-spin Kohn-Sham matrix.")
        .def(
            "alpha_to_numpy",
            [](const CAOKohnShamMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getPointerToAlphaValues(), {self.getNumberOfRows(), self.getNumberOfColumns()});
            },
            "Converts alpha AOKohnShamMatrix to numpy array.")
        .def(
            "beta_to_numpy",
            [](const CAOKohnShamMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getPointerToBetaValues(), {self.getNumberOfRows(), self.getNumberOfColumns()});
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
                return vlx_general::pointer_to_numpy(self.getCoordinatesX(), std::vector<int64_t>{self.getNumberOfGridPoints()});
            },
            "Gets X coordinates of grid as numpy array.")
        .def(
            "y_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesY(), std::vector<int64_t>{self.getNumberOfGridPoints()});
            },
            "Gets Y coordinates of grid as numpy array.")
        .def(
            "z_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), std::vector<int64_t>{self.getNumberOfGridPoints()});
            },
            "Gets Z coordinates of grid as numpy array.")
        .def(
            "w_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getWeights(), std::vector<int64_t>{self.getNumberOfGridPoints()});
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
        .def(py::init(&CXCIntegrator_create), "comm"_a = py::none())
        .def("integrate_vxc_fock",
             &CXCIntegrator::integrateVxcFock,
             "Integrates 1st-order exchange-correlation contribution to Kohn-Sham matrix.",
             "molecule"_a,
             "basis"_a,
             "density_matrix"_a,
             "molecular_grid"_a,
             "flag"_a)
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

    // XCComponent class
    PyClass<CXCComponent>(m, "XCComponent")
        .def(py::init<const std::string&, const double>(), "label"_a, "coeff"_a)
        .def(py::init<const CXCComponent&>())
        .def("get_scaling_factor", &CXCComponent::getScalingFactor, "Gets scaling factor of XC functional component.")
        .def("get_label", &CXCComponent::getLabel, "Gets name of XC functional component.")
        .def(py::self == py::self);

    // exposing functions

    m.def("available_functionals", &vxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
