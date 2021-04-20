//
//                           VELOXCHEM 1.0-RC
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

#include "ExportVisualization.hpp"

#include <mpi.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CubicGrid.hpp"
#include "ExportGeneral.hpp"
#include "VisualizationDriver.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_visualization {  // vlx_visualization namespace

// Helper function for converting cubic grid values to 3d numpy array

static py::array_t<double>
CCubicGrid_values_to_numpy(const CCubicGrid& self)
{
    auto nx = self.numPointsX();

    auto ny = self.numPointsY();

    auto nz = self.numPointsZ();

    return vlx_general::pointer_to_numpy(self.values(), {nx, ny, nz});
}

// Exports classes/functions in src/visualization to python

void
export_visualization(py::module& m)
{
    // CCubicGrid class

    PyClass<CCubicGrid>(m, "CubicGrid")
        .def(py::init<>())
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const std::vector<int32_t>>())
        .def("x_origin", &CCubicGrid::originX)
        .def("y_origin", &CCubicGrid::originY)
        .def("z_origin", &CCubicGrid::originZ)
        .def("x_step_size", &CCubicGrid::stepSizeX)
        .def("y_step_size", &CCubicGrid::stepSizeY)
        .def("z_step_size", &CCubicGrid::stepSizeZ)
        .def("x_num_points", &CCubicGrid::numPointsX)
        .def("y_num_points", &CCubicGrid::numPointsY)
        .def("z_num_points", &CCubicGrid::numPointsZ)
        .def("set_values", &CCubicGrid::setValues)
        .def("values_to_numpy", &CCubicGrid_values_to_numpy);

    // CVisualizationDriver class

    PyClass<CVisualizationDriver>(m, "VisualizationDriver")
        .def(py::init(&vlx_general::create<CVisualizationDriver>), "comm"_a = py::none())
        .def("get_rank", &CVisualizationDriver::getRank)
        .def(
            "compute",
            vlx_general::
                overload_cast_<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CMolecularOrbitals&, const int32_t, const std::string&>()(
                    &CVisualizationDriver::compute, py::const_),
            "Compute molecular orbital values on cubic grid",
            "grid"_a,
            "molecule"_a,
            "basis"_a,
            "mos"_a,
            "mo_idx"_a,
            "mo_spin"_a)
        .def("compute",
             vlx_general::
                 overload_cast_<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CAODensityMatrix&, const int32_t, const std::string&>()(
                     &CVisualizationDriver::compute, py::const_),
             "Compute density values on cubic grid",
             "grid"_a,
             "molecule"_a,
             "basis"_a,
             "ao_density"_a,
             "idx"_a,
             "spin"_a)
        .def("get_mo", &CVisualizationDriver::getMO)
        .def("get_density", &CVisualizationDriver::getDensity)
        .def("get_two_particle_density", &CVisualizationDriver::getTwoParticleDensity);
}

}  // namespace vlx_visualization
