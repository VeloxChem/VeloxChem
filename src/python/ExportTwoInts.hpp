//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportTwoInts_hpp
#define ExportTwoInts_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_twoints {  // vlx_twoints namespace

/**
 Exports classes/functions in src/twoints to python.
 */
void export_twoints(py::module& m);

}  // namespace vlx_twoints

#endif /* ExportTwoInts_hpp */
