//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportVisualization_hpp
#define ExportVisualization_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_visualization { // vlx_visualization namespace

/**
 Exports classes/functions in src/visualization to python.
 */
void export_visualization(py::module& m);

} // vlx_visualization namespace

#endif /* ExportVisualization_hpp */
