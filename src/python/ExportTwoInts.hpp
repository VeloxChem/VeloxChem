//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

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
