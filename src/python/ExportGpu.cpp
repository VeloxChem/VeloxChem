//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#ifdef ENABLE_GPU
#include "DeviceProp.hpp"
#endif

#include "ExportGpu.hpp"

namespace bp = boost::python;

namespace bp_gpu { // bp_gpu namespace

// Exports classes/functions in src/gpu to python

void export_gpu()
{
    #ifdef ENABLE_GPU

    bp::def("get_device_properties", &gpu::getDeviceProperties);

    #endif
}

} // bp_gpu namespace
