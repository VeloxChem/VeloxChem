//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>

#ifdef ENABLE_GPU
#include "DeviceProp.hpp"
#endif

#include "ExportGpu.hpp"

namespace py = pybind11;

namespace vlx_gpu {  // vlx_gpu namespace

// Exports classes/functions in src/gpu to python

void
export_gpu(py::module& m)
{
#ifdef ENABLE_GPU

    m.def("get_device_properties", &gpu::getDeviceProperties);

#endif
}

}  // namespace vlx_gpu
