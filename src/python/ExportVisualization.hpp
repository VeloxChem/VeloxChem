//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportVisualization_hpp
#define ExportVisualization_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_visualization {  // vlx_visualization namespace

/**
 Exports classes/functions in src/visualization to python.
 */
void export_visualization(py::module& m);

}  // namespace vlx_visualization

#endif /* ExportVisualization_hpp */
