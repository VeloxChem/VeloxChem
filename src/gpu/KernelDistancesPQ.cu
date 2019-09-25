//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelDistancesPQ.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelDistancesPQ(      double* pqDistancesData,
                                  const size_t  pitchOfDistancesData,
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
        // set up pointers to P center coordinates

        const double* rpx = (double*)((char*)braGtoPairsData + 4 * pitchOfBraGtoPairsData);

        const double* rpy = (double*)((char*)braGtoPairsData + 5 * pitchOfBraGtoPairsData);

        const double* rpz = (double*)((char*)braGtoPairsData + 6 * pitchOfBraGtoPairsData);

        // set up pointers to Q center coordinates

        const double* rqx = (double*)((char*)ketGtoPairsData + 4 * pitchOfKetGtoPairsData);

        const double* rqy = (double*)((char*)ketGtoPairsData + 5 * pitchOfKetGtoPairsData);

        const double* rqz = (double*)((char*)ketGtoPairsData + 6 * pitchOfKetGtoPairsData);

        // compute R(PQ) = P - Q distances

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to R(PQ) distances

            int32_t ioff = 3 * (i - startPositionOfBraPair);

            double* pqx = (double*)((char*)pqDistancesData + ioff * pitchOfDistancesData);

            double* pqy = (double*)((char*)pqDistancesData + (ioff + 1) * pitchOfDistancesData);

            double* pqz = (double*)((char*)pqDistancesData + (ioff + 2) * pitchOfDistancesData);

            // compute R(PQ) distances

            pqx[tid] = rpx[i] - rqx[tid];

            pqy[tid] = rpy[i] - rqy[tid];

            pqz[tid] = rpz[i] - rqz[tid];
        }
    }
}

void
launchKernelForDistancesPQ(      double* pqDistancesData,
                           const size_t  pitchOfDistancesData,
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
    gpu::kernelDistancesPQ<<<gridSize, blockSize>>>(pqDistancesData, pitchOfDistancesData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    startPositionOfBraPair, endPositionOfBraPair,
                                                    nKetPrimPairs);
}

}  // namespace gpu

