//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportOneInts_hpp
#define ExportOneInts_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_oneints {  // vlx_oneints namespace

/**
 Exports classes/functions in src/oneints to python.
 */
void export_oneints(py::module& m);

}  // namespace vlx_oneints

#endif /* ExportOneInts_hpp */
