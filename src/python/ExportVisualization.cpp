//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>

#include "VisualizationDriver.hpp"
#include "ExportVisualization.hpp"

namespace py = pybind11;

namespace bp_visualization { // bp_visualization namespace

// Exports classes/functions in src/visualization to python

void export_visualization(py::module& m)
{
    // CVisualizationDriver class

    py::class_< CVisualizationDriver, std::shared_ptr<CVisualizationDriver> >
        (
            m, "VisualizationDriver"
        )
        .def(py::init<>())
        .def("compute",
             (double (CVisualizationDriver::*)(const CMolecule&,
                                               const CMolecularBasis&,
                                               const CMolecularOrbitals&,
                                               const int32_t,
                                               const std::string&,
                                               const double,
                                               const double,
                                               const double) const)
             &CVisualizationDriver::compute)
        .def("compute",
             (double (CVisualizationDriver::*)(const CMolecule&,
                                               const CMolecularBasis&,
                                               const CAODensityMatrix&,
                                               const int32_t,
                                               const std::string&,
                                               const double,
                                               const double,
                                               const double) const)
             &CVisualizationDriver::compute)
    ;
}

} // bp_visualization namespace
