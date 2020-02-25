//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
