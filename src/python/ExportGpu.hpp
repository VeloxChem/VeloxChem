//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportGpu_hpp
#define ExportGpu_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_gpu { // vlx_gpu namespace

/**
 Exports classes/functions in src/gpu to python.
 */
void export_gpu(py::module& m);

} // vlx_gpu namespace

#endif /* ExportGpu_hpp */
