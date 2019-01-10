//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportOrbData_hpp
#define ExportOrbData_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_orbdata { // vlx_orbdata namespace

/**
 Exports classes/functions in src/orbdata to python.
 */
void export_orbdata(py::module& m);

} // vlx_orbdata namespace

#endif /* ExportOrbData_hpp */
