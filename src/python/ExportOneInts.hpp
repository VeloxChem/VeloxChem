//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

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
