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

        int32_t devcnt = 0;

        cudaGetDeviceCount(&devcnt);

        oStream << fmt::info;

        for (int32_t i = 0; i < devcnt; i++) {

            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            oStream << "GPU device ID: " << std::to_string(i) << fmt::end;

            oStream << "  Device name:             ";

            oStream << prop.name << fmt::end;

            oStream << "  Compute capability:      ";

            oStream << std::to_string(prop.major) << "." << std::to_string(prop.minor) << fmt::end;

            oStream << "  Multiprocessor count:    ";
            
            oStream << std::to_string(prop.multiProcessorCount) << fmt::end;

            oStream << "  Max clock rate:          ";
            
            oStream << fstr::to_string(prop.clockRate * 1.0e-6, 2) << " GHz" << fmt::end;

            double glbmem = (double)prop.totalGlobalMem / std::pow(1024, 3);

            oStream << "  Global memory:           ";
            
            oStream << fstr::to_string(glbmem, 0) << " GB" << fmt::end;

            oStream << "  Peak memory bandwidth:   ";

            double bandwidth = 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e+6;
            
            oStream << fstr::to_string(bandwidth, 0) << " GB/s" << fmt::end;

        }

        oStream << fmt::blank;

    }

    #endif

}
