//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "CubicGrid.hpp"
#include "VisualizationDriver.hpp"
#include "ExportGeneral.hpp"
#include "ExportVisualization.hpp"

namespace py = pybind11;

namespace vlx_visualization { // vlx_visualization namespace

// Helper functions for CVisualizationDriver::compute

static py::array_t<double>
CVisualizationDriver_compute_mo(      CVisualizationDriver& self,
                                const CMolecule&            molecule,
                                const CMolecularBasis&      basis,
                                const CMolecularOrbitals&   mo,
                                const int32_t               moidx,
                                const std::string&          mospin,
                                const CCubicGrid&           grid)
{
    auto psi_data = self.compute(molecule, basis, mo, moidx, mospin, grid);

    auto nx = grid.numPointsX();

    auto ny = grid.numPointsY();

    auto nz = grid.numPointsZ();

    return vlx_general::pointer_to_numpy(psi_data.data(), {nx, ny, nz});
}

static py::array_t<double>
CVisualizationDriver_compute_density(      CVisualizationDriver& self,
                                     const CMolecule&            molecule,
                                     const CMolecularBasis&      basis,
                                     const CAODensityMatrix&     density,
                                     const int32_t               denidx,
                                     const std::string&          denspin,
                                     const CCubicGrid&           grid)
{
    auto psi_data = self.compute(molecule, basis, density, denidx, denspin, grid);

    auto nx = grid.numPointsX();

    auto ny = grid.numPointsY();

    auto nz = grid.numPointsZ();

    return vlx_general::pointer_to_numpy(psi_data.data(), {nx, ny, nz});
}

// Exports classes/functions in src/visualization to python

void export_visualization(py::module& m)
{
    // CCubicGrid class

    py::class_< CCubicGrid, std::shared_ptr<CCubicGrid> >
        (
            m, "CubicGrid"
        )
        .def(py::init<>())
        .def(py::init<const std::vector<double>&,
                      const std::vector<double>&,
                      const std::vector<int32_t>>())
        .def("x_origin", &CCubicGrid::originX)
        .def("y_origin", &CCubicGrid::originY)
        .def("z_origin", &CCubicGrid::originZ)
        .def("x_step_size", &CCubicGrid::stepSizeX)
        .def("y_step_size", &CCubicGrid::stepSizeY)
        .def("z_step_size", &CCubicGrid::stepSizeZ)
        .def("x_num_points", &CCubicGrid::numPointsX)
        .def("y_num_points", &CCubicGrid::numPointsY)
        .def("z_num_points", &CCubicGrid::numPointsZ)
    ;

    // CVisualizationDriver class

    py::class_< CVisualizationDriver, std::shared_ptr<CVisualizationDriver> >
        (
            m, "VisualizationDriver"
        )
        .def(py::init<>())
        .def("compute", &CVisualizationDriver_compute_mo)
        .def("compute", &CVisualizationDriver_compute_density)
    ;
}

} // vlx_visualization namespace
