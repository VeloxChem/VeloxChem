//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>

#include "CudaDevices.hpp"

#include "ExportGpu.hpp"

namespace py = pybind11;

namespace vlx_gpu {  // vlx_gpu namespace

// Helper function for printing CCudaDevices
    
static std::string
CCudaDevices_str(const CCudaDevices& self)
{
    return self.getString();
}
    
// Exports classes/functions in src/gpu to python

void
export_gpu(py::module& m)
{
    // CCudaDevices class
    
    py::class_<CCudaDevices, std::shared_ptr<CCudaDevices>>(m, "CudaDevices")
        .def(py::init<>())
        .def("get_number_devices", &CCudaDevices::getNumberOfDevices)
        .def("__str__", &CCudaDevices_str);
}

}  // namespace vlx_gpu
