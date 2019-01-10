//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportMolData_hpp
#define ExportMolData_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_moldata { // vlx_moldata namespace

/**
 Exports classes/functions in src/moldata to python.
 */
void export_moldata(py::module& m);

} // vlx_moldata namespace

#endif /* ExportMolData_hpp */
