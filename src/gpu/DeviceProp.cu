//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cmath>
#include <string>
#include <sstream>

#include "DeviceProp.hpp"
#include "StringFormat.hpp"

namespace gpu { // gpu namespace

    std::string getDeviceProperties()
    {
        int32_t devcnt = 0;

        cudaGetDeviceCount(&devcnt);

        std::stringstream ss;

        for (int32_t i = 0; i < devcnt; i++)
        {
            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            ss << "GPU device ID: " << std::to_string(i) << "\n";

            ss << "  Device name:             ";

            ss << prop.name << "\n";

            ss << "  Compute capability:      ";

            ss << std::to_string(prop.major) << "." << std::to_string(prop.minor) << "\n";

            ss << "  Multiprocessor count:    ";
            
            ss << std::to_string(prop.multiProcessorCount) << "\n";

            ss << "  Max clock rate:          ";
            
            ss << fstr::to_string(prop.clockRate * 1.0e-6, 2) << " GHz" << "\n";

            double glbmem = (double)prop.totalGlobalMem / std::pow(1024, 3);

            ss << "  Global memory:           ";
            
            ss << fstr::to_string(glbmem, 0) << " GB" << "\n";

            ss << "  Peak memory bandwidth:   ";

            double bandwidth = 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e+6;
            
            ss << fstr::to_string(bandwidth, 0) << " GB/s" << "\n\n";
        }

        return ss.str();
    }

}
