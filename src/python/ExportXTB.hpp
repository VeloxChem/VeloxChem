//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportXTB_hpp
#define ExportXTB_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_xtb {  // vlx_xtb namespace

/**
 Exports classes/functions in src/xtb to python.
 */
void export_xtb(py::module& m);

}  // namespace vlx_xtb

#endif /* ExportXTB_hpp */
