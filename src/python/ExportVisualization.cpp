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
        .def("x_origin", &CCubicGrid::originX, "Gets X coordinate of the origin.")
        .def("y_origin", &CCubicGrid::originY, "Gets Y coordinate of the origin.")
        .def("z_origin", &CCubicGrid::originZ, "Gets Z coordinate of the origin.")
        .def("x_step_size", &CCubicGrid::stepSizeX, "Gets step size in X direction.")
        .def("y_step_size", &CCubicGrid::stepSizeY, "Gets step size in Y direction.")
        .def("z_step_size", &CCubicGrid::stepSizeZ, "Gets step size in Z direction.")
        .def("x_num_points", &CCubicGrid::numPointsX, "Gets number of points in X direction.")
        .def("y_num_points", &CCubicGrid::numPointsY, "Gets number of points in Y direction.")
        .def("z_num_points", &CCubicGrid::numPointsZ, "Gets number of points in Z direction.")
        .def("set_values", &CCubicGrid::setValues, "Sets the cubic grid values.", "vals"_a)
        .def("values_to_numpy", &CCubicGrid_values_to_numpy, "Convertis cubic grid values to 3D numpy array.");

    // CVisualizationDriver class

    PyClass<CVisualizationDriver>(m, "VisualizationDriver")
        .def(py::init(&vlx_general::create<CVisualizationDriver>), "comm"_a = py::none())
        .def("get_rank", &CVisualizationDriver::getRank, "Gets rank of the MPI process.")
        .def(
            "compute",
            vlx_general::
                overload_cast_<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CMolecularOrbitals&, const int32_t, const std::string&>()(
                    &CVisualizationDriver::compute, py::const_),
            "Computes molecular orbital values at cubic grid points.",
            "grid"_a,
            "molecule"_a,
            "basis"_a,
            "molorb"_a,
            "moidx"_a,
            "mospin"_a)
        .def("compute",
             vlx_general::
                 overload_cast_<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CAODensityMatrix&, const int32_t, const std::string&>()(
                     &CVisualizationDriver::compute, py::const_),
             "Computes density values at cubic grid points.",
             "grid"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "denidx"_a,
             "denspin"_a)
        .def("get_mo",
             &CVisualizationDriver::getMO,
             "Computes molecular orbital at given coordinates.",
             "coords"_a,
             "molecule"_a,
             "basis"_a,
             "mo"_a,
             "moidx"_a,
             "mospin"_a)
        .def("get_density",
             &CVisualizationDriver::getDensity,
             "Computes densities at given coordinates.",
             "coords"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "denidx"_a,
             "denspin"_a)
        .def("get_two_particle_density",
             &CVisualizationDriver::getTwoParticleDensity,
             "Computes two-particle density Gamma(x1,x2;x1,x2) at given coordinates.",
             "coords_1"_a,
             "coords_2"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "denidx"_a,
             "spin_1"_a,
             "spin_2"_a);
}

}  // namespace vlx_visualization
