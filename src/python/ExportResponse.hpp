//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportResponse_hpp
#define ExportResponse_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_response {  // vlx_response namespace

/**
 Exports classes/functions in src/response to python.
*/
void export_response(py::module& m);

}  // namespace vlx_response

#endif /* ExportResponse_hpp */
