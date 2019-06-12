//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportMolData_hpp
#define ExportMolData_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_moldata {  // vlx_moldata namespace

/**
 Exports classes/functions in src/moldata to python.
 */
void export_moldata(py::module& m);

}  // namespace vlx_moldata

#endif /* ExportMolData_hpp */
