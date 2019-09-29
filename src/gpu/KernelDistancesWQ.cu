//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelDistancesWQ.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelDistancesWQ(      double* wqDistancesData,
                                  const size_t  pitchOfDistancesData,
                                  const double* wCoordinatesData,
                                  const size_t  pitchOfCoordinatesData,
                                  const double* ketGtoPairsData,
                                  const size_t  pitchOfKetGtoPairsData,
                                  const int32_t startPositionOfBraPair,
                                  const int32_t endPositionOfBraPair,
                                  const int32_t nKetPrimPairs)
{
    int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < nKetPrimPairs)
    {
        // set up pointers to Q center coordinates

        const double* rqx = (double*)((char*)ketGtoPairsData + 4 * pitchOfKetGtoPairsData);

        const double* rqy = (double*)((char*)ketGtoPairsData + 5 * pitchOfKetGtoPairsData);

        const double* rqz = (double*)((char*)ketGtoPairsData + 6 * pitchOfKetGtoPairsData);

        // compute distances R(WQ) = W - Q

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to coordinates W

            int32_t ioff = 3 * (i - startPositionOfBraPair);

            const double* wx = (double*)((char*)wCoordinatesData + ioff * pitchOfCoordinatesData);

            const double* wy = (double*)((char*)wCoordinatesData + (ioff + 1) * pitchOfCoordinatesData);

            const double* wz = (double*)((char*)wCoordinatesData + (ioff + 2) * pitchOfCoordinatesData);

            // set up pointers to distances R(WQ) = W - Q

            double* wqx = (double*)((char*)wqDistancesData + ioff * pitchOfDistancesData);

            double* wqy = (double*)((char*)wqDistancesData + (ioff + 1) * pitchOfDistancesData);

            double* wqz = (double*)((char*)wqDistancesData + (ioff + 2) * pitchOfDistancesData);

            // compute coordinates R(WQ) = W - Q

            wqx[tid] = wx[tid] - rqx[tid];

            wqy[tid] = wy[tid] - rqy[tid];

            wqz[tid] = wz[tid] - rqz[tid];
        }
    }
}

void launchKernelForDistancesWQ(      double* wqDistancesData,
                                const size_t  pitchOfDistancesData,
                                const double* wCoordinatesData,
                                const size_t  pitchOfCoordinatesData,
                                const double* ketGtoPairsData,
                                const size_t  pitchOfKetGtoPairsData,
                                const int32_t startPositionOfBraPair,
                                const int32_t endPositionOfBraPair,
                                const int32_t nKetPrimPairs,
                                const int32_t gridSize,
                                const int32_t blockSize)
{
    gpu::kernelDistancesWQ<<<gridSize, blockSize>>>(wqDistancesData, pitchOfDistancesData,
                                                    wCoordinatesData, pitchOfCoordinatesData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    startPositionOfBraPair, endPositionOfBraPair,
                                                    nKetPrimPairs);
}

}  // namespace gpu

