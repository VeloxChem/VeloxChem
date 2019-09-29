//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelDistancesWP.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelDistancesWP(      double* wpDistancesData,
                                  const size_t  pitchOfDistancesData,
                                  const double* wCoordinatesData,
                                  const size_t  pitchOfCoordinatesData,
                                  const double* braGtoPairsData,
                                  const size_t  pitchOfBraGtoPairsData,
                                  const int32_t startPositionOfBraPair,
                                  const int32_t endPositionOfBraPair,
                                  const int32_t nKetPrimPairs)
{
    int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < nKetPrimPairs)
    {
        // set up pointers to P center coordinates

        const double* rpx = (double*)((char*)braGtoPairsData + 4 * pitchOfBraGtoPairsData);

        const double* rpy = (double*)((char*)braGtoPairsData + 5 * pitchOfBraGtoPairsData);

        const double* rpz = (double*)((char*)braGtoPairsData + 6 * pitchOfBraGtoPairsData);

        // compute distances R(WP) = W - P

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to coordinates W

            int32_t ioff = 3 * (i - startPositionOfBraPair);

            const double* wx = (double*)((char*)wCoordinatesData + ioff * pitchOfCoordinatesData);

            const double* wy = (double*)((char*)wCoordinatesData + (ioff + 1) * pitchOfCoordinatesData);

            const double* wz = (double*)((char*)wCoordinatesData + (ioff + 2) * pitchOfCoordinatesData);

            // set up pointers to distances R(WP) = W - P

            double* wpx = (double*)((char*)wpDistancesData + ioff * pitchOfDistancesData);

            double* wpy = (double*)((char*)wpDistancesData + (ioff + 1) * pitchOfDistancesData);

            double* wpz = (double*)((char*)wpDistancesData + (ioff + 2) * pitchOfDistancesData);

            // compute coordinates R(WP) = W - P

            wpx[tid] = wx[tid] - rpx[i];

            wpy[tid] = wy[tid] - rpy[i];

            wpz[tid] = wz[tid] - rpz[i];
        }
    }
}

void launchKernelForDistancesWP(      double* wpDistancesData,
                                const size_t  pitchOfDistancesData,
                                const double* wCoordinatesData,
                                const size_t  pitchOfCoordinatesData,
                                const double* braGtoPairsData,
                                const size_t  pitchOfBraGtoPairsData,
                                const int32_t startPositionOfBraPair,
                                const int32_t endPositionOfBraPair,
                                const int32_t nKetPrimPairs,
                                const int32_t gridSize,
                                const int32_t blockSize)
{
    gpu::kernelDistancesWP<<<gridSize, blockSize>>>(wpDistancesData, pitchOfDistancesData,
                                                    wCoordinatesData, pitchOfCoordinatesData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    startPositionOfBraPair, endPositionOfBraPair,
                                                    nKetPrimPairs);
}

}  // namespace gpu

