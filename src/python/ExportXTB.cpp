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
#include "XTBDriver.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_xtb {  // vlx_xtb namespace

// Exports classes/functions in src/xtb to python

void
export_xtb(py::module& m)
{
    // CXTBDriver class

    PyClass<CXTBDriver>(m, "XTBDriver")
        .def(py::init(&vlx_general::create<CXTBDriver>), "comm"_a = py::none())
        .def_static("is_available", &CXTBDriver::isAvailable, "Checks if XTB package is available.")
        .def("is_master_node", &CXTBDriver::isMasterNode, "Checks if XTB driver is running on master node.")
        .def("set_max_iter", &CXTBDriver::setMaxIterations, "Sets maximum number of SCF iterations.", "maxIterations"_a)
        .def("set_elec_temp", &CXTBDriver::setElectronicTemp, "Sets electronic temperature for electron smearing.", "electronicTemp"_a)
        .def("set_method", &CXTBDriver::setMethod, "Sets XTB method.", "method"_a)
        .def("set_output_filename", &CXTBDriver::setOutputFilename, "Sets output filename.", "filename"_a)
        .def("get_output", &CXTBDriver::getOutput, "Gets XTB output as a vector of strings.")
        .def("get_output_filename", &CXTBDriver::getOutputFilename, "Gets XTB output filename.")
        .def("get_energy", &CXTBDriver::getEnergy, "Gets energy of molecular system.")
        .def(
            "get_gradient",
            [](const CXTBDriver& self) -> py::array_t<double> {
                auto grad = self.getGradient();
                return vlx_general::pointer_to_numpy(grad.data(), grad.size() / 3, 3);
            },
            "Gets molecular gradient as numpy array of shape (natoms, 3).")
        .def(
            "get_dipole",
            [](const CXTBDriver& self) -> py::array_t<double> {
                auto dipole = self.getDipole();
                return vlx_general::pointer_to_numpy(dipole.data(), dipole.size());
            },
            "Gets molecular dipole moment as numpy array.")
        // prefixed by an underscore because it will be decorated in xtbdriver.py
        .def("_compute", &CXTBDriver::compute);
}

}  // namespace vlx_xtb
