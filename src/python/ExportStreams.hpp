//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportStreams_hpp
#define ExportStreams_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_streams { // vlx_streams namespace

/**
 Exports classes/functions in src/streams to python.
 */
void export_streams(py::module& m);

} // vlx_streams namespace

#endif /* ExportStreams_hpp */
