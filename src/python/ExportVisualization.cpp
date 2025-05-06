//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ExportVisualization.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <vector>

#include "AtomicRadii.hpp"
#include "CubicGrid.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "VisualizationDriver.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_visualization {  // vlx_visualization namespace

static auto
CVisualizationDriver_compute_local_grid(CVisualizationDriver&      self,
                             CCubicGrid&                grid,
                             const CMolecule&           molecule,
                             const CMolecularBasis&     basis,
                             const py::array_t<double>& mocoefs,
                             const int                  moidx) -> void
{
    auto nao = mocoefs.shape(0);
    auto nmo = mocoefs.shape(1);

    self.compute_local_grid(grid, molecule, basis, nao, nmo, mocoefs.data(), moidx);
}

static auto
CVisualizationDriver_get_mo(CVisualizationDriver&                   self,
                            const std::vector<std::vector<double>>& coords,
                            const CMolecule&                        molecule,
                            const CMolecularBasis&                  basis,
                            const py::array_t<double>&              mocoefs,
                            const int                               moidx) -> std::vector<double>
{
    auto nao = mocoefs.shape(0);
    auto nmo = mocoefs.shape(1);

    return self.getMO(coords, molecule, basis, nao, nmo, mocoefs.data(), moidx);
}

// Exports classes/functions in src/visualization to python

void
export_visualization(py::module& m)
{
    // CCubicGrid class

    PyClass<CCubicGrid>(m, "CubicGrid")
        .def(py::init<>())
        .def(py::init<const std::array<double, 3>&, const std::array<double, 3>&, const std::array<int, 3>>())
        .def("get_origin", &CCubicGrid::getOrigin, "Gets coordinate of the origin.")
        .def("get_step_size", &CCubicGrid::getStepSize, "Gets step size in X, Y and Z direction.")
        .def("get_num_points", &CCubicGrid::getNumPoints, "Gets number of points in X, Y and Z direction.")
        .def("set_values", &CCubicGrid::setValues, "Sets the cubic grid values.", "vals"_a)
        .def(
            "values_to_numpy",
            [](const CCubicGrid& self) -> py::array_t<double> {
              auto num_points = self.getNumPoints();
              return vlx_general::pointer_to_numpy(self.values(), {num_points[0], num_points[1], num_points[2]});
            },
            "Convertis cubic grid values to 3D numpy array.");

    // CVisualizationDriver class
    // Note: VisualizationDriver is prefixed by an underscore and will be used in visualizationdriver.py

    PyClass<CVisualizationDriver>(m, "_VisualizationDriver")
        .def(py::init<>())
        .def("_create_local_cubic_grid", &CVisualizationDriver::create_local_cubic_grid, "Creates MPI-local cubic grid.")
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
        .def("_compute_local_grid",
             &CVisualizationDriver_compute_local_grid,
             "Computes molecular orbital values at MPI-local cubic grid points.",
             "grid"_a,
             "molecule"_a,
             "basis"_a,
             "mocoefs"_a,
             "moidx"_a)
        .def("_compute_local_grid",
             py::overload_cast<CCubicGrid&, const CMolecule&, const CMolecularBasis&, const CAODensityMatrix&, const int, const std::string&>(
                 &CVisualizationDriver::compute_local_grid, py::const_),
             "Computes density values at MPI-local cubic grid points.",
             "grid"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "denidx"_a,
             "denspin"_a)
        .def("get_mo",
             &CVisualizationDriver_get_mo,
             "Computes molecular orbital at given coordinates.",
             "coords"_a,
             "molecule"_a,
             "basis"_a,
             "mocoefs"_a,
             "moidx"_a)
        .def("get_density",
             &CVisualizationDriver::getDensity,
             "Computes densities at given coordinates.",
             "coords"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "denspin"_a)
        .def("get_two_particle_density",
             &CVisualizationDriver::getTwoParticleDensity,
             "Computes two-particle density Gamma(x1,x2;x1,x2) at given coordinates.",
             "coords_1"_a,
             "coords_2"_a,
             "molecule"_a,
             "basis"_a,
             "density"_a,
             "spin_1"_a,
             "spin_2"_a);
}

}  // namespace vlx_visualization
