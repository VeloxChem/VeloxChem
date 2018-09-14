//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cmath>
#include <string>

#include "DeviceProp.hpp"
#include "StringFormat.hpp"

namespace gpu { // gpu namespace

    #ifdef ENABLE_GPU

    void get_device_prop(COutputStream& oStream)
    {

        int32_t nDevices = 0;

        cudaGetDeviceCount(&nDevices);

        oStream << fmt::info << "Total Number of GPU Devices: ";

        oStream << std::to_string(nDevices) << fmt::end;

        for (int32_t i = 0; i < nDevices; i++) {

            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            oStream << "Device ID: " << std::to_string(i) << fmt::end;

            oStream << "  Device name:             ";

            oStream << prop.name << fmt::end;

            oStream << "  Compute Capability:      ";

            oStream << std::to_string(prop.major) << "." << std::to_string(prop.minor) << fmt::end;

            oStream << "  Multiprocessor Count:    ";
            
            oStream << std::to_string(prop.multiProcessorCount) << fmt::end;

            oStream << "  Max Clock Rate:          ";
            
            oStream << fstr::to_string(prop.clockRate * 1.0e-6, 2) << " GHz" << fmt::end;

            double globalMem = (double)prop.totalGlobalMem / std::pow(1024, 3);

            oStream << "  Global Memory:           ";
            
            oStream << fstr::to_string(globalMem, 0) << " GB" << fmt::end;

            oStream << "  Peak Memory Bandwidth:   ";

            double memBandwidth = 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e+6;
            
            oStream << fstr::to_string(memBandwidth, 0) << " GB/s" << fmt::end;

        }

    }

    #endif

}
