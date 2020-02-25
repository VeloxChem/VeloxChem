//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportOrbData_hpp
#define ExportOrbData_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_orbdata {  // vlx_orbdata namespace

/**
 Exports classes/functions in src/orbdata to python.
 */
void export_orbdata(py::module& m);

}  // namespace vlx_orbdata

#endif /* ExportOrbData_hpp */
