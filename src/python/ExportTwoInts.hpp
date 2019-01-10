//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportTwoInts_hpp
#define ExportTwoInts_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_twoints { // vlx_twoints namespace

/**
 Exports classes/functions in src/twoints to python.
 */
void export_twoints(py::module& m);

} // vlx_twoints namespace

#endif /* ExportTwoInts_hpp */
