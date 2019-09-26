//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelCoordinatesW.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelCoordinatesW(      double* wCoordinatesData,
                                   const size_t  pitchOfCoordinatesData,
                                   const double* osFactorsData,
                                   const size_t  pitchOfFactorsData,
                                   const double* braGtoPairsData,
                                   const size_t  pitchOfBraGtoPairsData,
                                   const double* ketGtoPairsData,
                                   const size_t  pitchOfKetGtoPairsData,
                                   const int32_t startPositionOfBraPair,
                                   const int32_t endPositionOfBraPair,
                                   const int32_t nKetPrimPairs)
{
    int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < nKetPrimPairs)
    {
        // set up pointers to primitive factors

        const double* bfxi = braGtoPairsData;

        const double* kfxi = ketGtoPairsData;

        // set up pointers to P center coordinates

        const double* rpx = (double*)((char*)braGtoPairsData + 4 * pitchOfBraGtoPairsData);

        const double* rpy = (double*)((char*)braGtoPairsData + 5 * pitchOfBraGtoPairsData);

        const double* rpz = (double*)((char*)braGtoPairsData + 6 * pitchOfBraGtoPairsData);

        // set up pointers to Q center coordinates

        const double* rqx = (double*)((char*)ketGtoPairsData + 4 * pitchOfKetGtoPairsData);

        const double* rqy = (double*)((char*)ketGtoPairsData + 5 * pitchOfKetGtoPairsData);

        const double* rqz = (double*)((char*)ketGtoPairsData + 6 * pitchOfKetGtoPairsData);

        // compute coordinates W

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to Obara-Saika factors

            int32_t joff = 4 * (i - startPositionOfBraPair);

            const double* fact = (double*)((char*)osFactorsData + joff * pitchOfFactorsData);

            // set up pointers to coordinates W

            int32_t ioff = 3 * (i - startPositionOfBraPair);

            double* wx = (double*)((char*)wCoordinatesData + ioff * pitchOfCoordinatesData);

            double* wy = (double*)((char*)wCoordinatesData + (ioff + 1) * pitchOfCoordinatesData);

            double* wz = (double*)((char*)wCoordinatesData + (ioff + 2) * pitchOfCoordinatesData);

            // compute coordinates W

            const double bfact = bfxi[i];

            wx[tid] = fact[tid] * (bfact * rpx[i] + kfxi[tid] * rqx[tid]);

            wy[tid] = fact[tid] * (bfact * rpy[i] + kfxi[tid] * rqy[tid]);

            wz[tid] = fact[tid] * (bfact * rpz[i] + kfxi[tid] * rqz[tid]);
        }
    }
}

void
launchKernelForCoordinatesW(      double* wCoordinatesData,
                            const size_t  pitchOfCoordinatesData,
                            const double* osFactorsData,
                            const size_t  pitchOfFactorsData,
                            const double* braGtoPairsData,
                            const size_t  pitchOfBraGtoPairsData,
                            const double* ketGtoPairsData,
                            const size_t  pitchOfKetGtoPairsData,
                            const int32_t startPositionOfBraPair,
                            const int32_t endPositionOfBraPair,
                            const int32_t nKetPrimPairs,
                            const int32_t gridSize,
                            const int32_t blockSize)
{
    gpu::kernelCoordinatesW<<<gridSize, blockSize>>>(wCoordinatesData, pitchOfCoordinatesData,
                                                     osFactorsData, pitchOfFactorsData,
                                                     braGtoPairsData, pitchOfBraGtoPairsData,
                                                     ketGtoPairsData, pitchOfKetGtoPairsData,
                                                     startPositionOfBraPair, endPositionOfBraPair,
                                                     nKetPrimPairs);
}

}  // namespace gpu

