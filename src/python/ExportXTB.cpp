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

#include "ExportXTB.hpp"

#include <mpi.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "ExportGeneral.hpp"
#include "Molecule.hpp"
#include "XtbDriver.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_xtb {  // vlx_xtb namespace

// Exports classes/functions in src/xtb to python

void
export_xtb(py::module& m)
{
    // CXtbDriver class
    // Note: XtbDriver is prefixed by an underscore and will be used in xtbdriver.py

    PyClass<CXtbDriver>(m, "_XtbDriver")
        .def(py::init(&vlx_general::create<CXtbDriver>), "comm"_a = py::none())
        .def_static("is_available", &CXtbDriver::isAvailable, "Checks if XTB package is available.")
        .def("set_max_iter", &CXtbDriver::setMaxIterations, "Sets maximum number of SCF iterations.", "maxIterations"_a)
        .def("set_elec_temp", &CXtbDriver::setElectronicTemp, "Sets electronic temperature for electron smearing.", "electronicTemp"_a)
        .def("set_method", &CXtbDriver::setMethod, "Sets XTB method.", "method"_a)
        .def("set_output_filename", &CXtbDriver::setOutputFilename, "Sets output filename.", "filename"_a)
        .def("mute", &CXtbDriver::mute, "Mutes output.")
        .def("unmute", &CXtbDriver::unmute, "Unmutes output.")
        .def("get_output", &CXtbDriver::getOutput, "Gets XTB output as a vector of strings.")
        .def("get_method", &CXtbDriver::getMethod, "Gets XTB method.")
        .def("get_output_filename", &CXtbDriver::getOutputFilename, "Gets XTB output filename.")
        .def("get_energy", &CXtbDriver::getEnergy, "Gets energy of molecular system.")
        .def(
            "get_gradient",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto grad = self.getGradient();
                return vlx_general::pointer_to_numpy(grad.data(), grad.size() / 3, 3);
            },
            "Gets molecular gradient as numpy array of shape (natoms, 3).")
        .def(
            "get_dipole",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto dipole = self.getDipole();
                return vlx_general::pointer_to_numpy(dipole.data(), dipole.size());
            },
            "Gets molecular dipole moment as numpy array.")
        .def(
            "get_partial_charges",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto partial_charges = self.getPartialCharges();
                return vlx_general::pointer_to_numpy(partial_charges.data(), partial_charges.size());
            },
            "Gets partial charges as numpy array.")
        .def(
            "get_bond_orders",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto n_atoms     = self.getNumberOfAtoms();
                auto bond_orders = self.getBondOrders();
                return vlx_general::pointer_to_numpy(bond_orders.data(), n_atoms, n_atoms);
            },
            "Gets bond orders as numpy array of shape (natoms, natoms).")
        .def(
            "get_orbital_energies",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto nao              = self.getNumberOfAOs();
                auto orbital_energies = self.getOrbitalEnergies();
                return vlx_general::pointer_to_numpy(orbital_energies.data(), nao);
            },
            "Gets orbital energies as numpy array.")
        .def(
            "get_orbital_occupations",
            [](const CXtbDriver& self) -> py::array_t<double> {
                auto nao                 = self.getNumberOfAOs();
                auto orbital_occupations = self.getOrbitalOccupations();
                return vlx_general::pointer_to_numpy(orbital_occupations.data(), nao);
            },
            "Gets orbital occupation numbers as numpy array.")
        .def("compute", &CXtbDriver::compute);
}

}  // namespace vlx_xtb
