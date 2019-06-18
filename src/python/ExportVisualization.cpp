//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CubicGrid.hpp"
#include "ExportGeneral.hpp"
#include "ExportVisualization.hpp"
#include "VisualizationDriver.hpp"

namespace py = pybind11;

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

    py::class_<CCubicGrid, std::shared_ptr<CCubicGrid>>(m, "CubicGrid")
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
        .def("values_to_numpy", &CCubicGrid_values_to_numpy);

    // CVisualizationDriver class

    py::class_<CVisualizationDriver, std::shared_ptr<CVisualizationDriver>>(m, "VisualizationDriver")
        .def(py::init<>())
        .def("compute",
             (void (CVisualizationDriver::*)(CCubicGrid&,
                                             const CMolecule&,
                                             const CMolecularBasis&,
                                             const CMolecularOrbitals&,
                                             const int32_t,
                                             const std::string&) const) &
                 CVisualizationDriver::compute)
        .def("compute",
             (void (CVisualizationDriver::*)(CCubicGrid&,
                                             const CMolecule&,
                                             const CMolecularBasis&,
                                             const CAODensityMatrix&,
                                             const int32_t,
                                             const std::string&) const) &
                 CVisualizationDriver::compute);
}

}  // namespace vlx_visualization
