//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportResponse_hpp
#define ExportResponse_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_response { // vlx_response namespace
    
/**
 Exports classes/functions in src/response to python.
*/
void export_response(py::module& m);
        
} // vlx_response namespace

#endif /* ExportResponse_hpp */
