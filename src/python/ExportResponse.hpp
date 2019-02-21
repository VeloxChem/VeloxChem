//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportResponse_hpp
#define ExportResponse_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_rsp { // vlx_rsp namespace
    
/**
 Exports classes/functions in src/response to python.
*/
void export_response(py::module& m);
        
} // vlx_rsp namespace

#endif /* ExportResponse_hpp */
