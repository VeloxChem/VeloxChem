//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cstdio>

#include "DeviceProp.hpp"

namespace gpu { // gpu namespace

    void get_device_prop()
    {
        #ifdef ENABLE_GPU

        int nDevices;

        cudaGetDeviceCount(&nDevices);

        printf("Total Number of Devices: %d\n", nDevices);

        for (int i = 0; i < nDevices; i++) {

            cudaDeviceProp prop;

            cudaGetDeviceProperties(&prop, i);

            printf("Device number: %d\n", i);

            printf("  Device name: %s\n", prop.name);

            printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);

            printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);

            printf("  Peak Memory Bandwidth (GB/s): %f\n",
                    2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
        }

        #endif
    }

}
