//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cmath>
#include <cstdio>

#include "DeviceProp.hpp"

namespace gpu { // gpu namespace

    #ifdef ENABLE_GPU

    void get_device_prop()
    {

        int nDevices = 0;

        cudaGetDeviceCount(&nDevices);

        printf("Total Number of Devices: %d\n", nDevices);

        for (int i = 0; i < nDevices; i++) {

            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            printf("Device ID: %d\n", i);

            printf("  Device name:             %s\n", prop.name);

            printf("  Compute Capability:      %d.%d\n", prop.major, prop.minor);

            printf("  Multiprocessor Count:    %d\n", prop.multiProcessorCount);

            printf("  Max Clock Rate:          %.2f GHz\n", prop.clockRate * 1.0e-6);

            printf("  Global Memory:           %.0f GB\n", (float)prop.totalGlobalMem/pow(1024,3));

            printf("  Peak Memory Bandwidth:   %.0f GB/s\n",
                    2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e+6);

        }

    }

    #endif

}
