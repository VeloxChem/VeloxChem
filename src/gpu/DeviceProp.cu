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
        std::string str("GPU Devices");

        std::stringstream ss;

        const int32_t width = 50;

        ss << str << "\n";

        ss << std::string(str.size() + 2, '=') << "\n\n";

        int32_t devcnt = 0;

        cudaGetDeviceCount(&devcnt);

        for (int32_t i = 0; i < devcnt; i++)
        {
            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            str.assign("GPU device ID: ");

            str.append(std::to_string(i));

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Device name:             ");

            str.append(prop.name);

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Compute capability:      ");

            str.append(std::to_string(prop.major));

            str.append(".");

            str.append(std::to_string(prop.minor));

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Multiprocessor count:    ");

            str.append(std::to_string(prop.multiProcessorCount));

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Max clock rate:          ");

            str.append(fstr::to_string(prop.clockRate * 1.0e-6, 2));

            str.append(" GHz");

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Global memory:           ");

            double glbmem = (double)prop.totalGlobalMem / std::pow(1024, 3);

            str.append(fstr::to_string(glbmem, 0));

            str.append(" GB");

            ss << fstr::format(str, width, fmt::left) << "\n";

            str.assign("  Peak memory bandwidth:   ");

            double bandwidth = 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e+6;

            str.append(fstr::to_string(bandwidth, 0));

            str.append(" GB/s");

            ss << fstr::format(str, width, fmt::left) << "\n";

            if (i < devcnt - 1) ss << "\n";
        }

        return ss.str();
    }

}
