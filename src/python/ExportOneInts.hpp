//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportOneInts_hpp
#define ExportOneInts_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_oneints { // vlx_oneints namespace

/**
 Exports classes/functions in src/oneints to python.
 */
void export_oneints(py::module& m);

} // vlx_oneints namespace

#endif /* ExportOneInts_hpp */
