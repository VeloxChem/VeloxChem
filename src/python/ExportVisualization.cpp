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

// Exports classes/functions in src/visualization to python

void
export_visualization(py::module& m)
{
    // CCubicGrid class

    PyClass<CCubicGrid>(m, "CubicGrid")
        .def(py::init<>())
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const std::vector<int32_t>>())
        .def("get_origin", &CCubicGrid::getOrigin, "Gets coordinate of the origin.")
        .def("get_step_size", &CCubicGrid::getStepSize, "Gets step size in X, Y and Z direction.")
        .def("get_num_points", &CCubicGrid::getNumPoints, "Gets number of points in X, Y and Z direction.")
        .def("set_values", &CCubicGrid::setValues, "Sets the cubic grid values.", "vals"_a)
        .def(
            "values_to_numpy",
            [](const CCubicGrid& self) -> py::array_t<double> { return vlx_general::pointer_to_numpy(self.values(), self.getNumPoints()); },
            "Convertis cubic grid values to 3D numpy array.");

    // CVisualizationDriver class

    PyClass<CVisualizationDriver>(m, "VisualizationDriver")
        .def(py::init(&vlx_general::create<CVisualizationDriver>), "comm"_a = py::none())
        .def("get_rank", &CVisualizationDriver::getRank, "Gets rank of the MPI process.")
        .def("get_atomic_orbital_info",
             &CVisualizationDriver::getAtomicOrbitalInformation,
             "Gets atomic orbital information.",
             "molecule"_a,
             "basis"_a)
        .def("map_atom_to_atomic_orbitals", &CVisualizationDriver::mapAtomToAtomicOrbitals, "Maps atom to atomic orbitals.", "molecule"_a, "basis"_a)
        .def("compute_atomic_orbital_for_grid",
             &CVisualizationDriver::computeAtomicOrbitalForGrid,
             "Computes atomic orbital (centered at origin) at cubic grid points.",
             "grid"_a,
             "basis"_a,
             "aoinfo"_a)
        .def("compute",
             py::overload_cast<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CMolecularOrbitals&, const int32_t, const std::string&>(
                 &CVisualizationDriver::compute, py::const_),
             "Computes molecular orbital values at cubic grid points.",
             "grid"_a,
             "molecule"_a,
             "basis"_a,
             "molorb"_a,
             "moidx"_a,
             "mospin"_a)
        .def("compute",
             py::overload_cast<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CAODensityMatrix&, const int32_t, const std::string&>(
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
