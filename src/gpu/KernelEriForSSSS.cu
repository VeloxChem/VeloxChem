//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelCoordinatesW.hpp"

#include <math.h>

#include  "CudaBoysFunction.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelEriForSSSS(      double* primBufferData,
                                 const size_t  pitchOfBufferData,
                                 const int32_t posIntegralInBuffer,
                                 const int32_t maxOrderOfIntegral,
                                 const double* osFactorsData,
                                 const size_t  pitchOfFactorsData,
                                 const double* pqDistancesData,
                                 const size_t  pitchOfDistancesPQData,
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
        // number of primitive integrals components

        int32_t ncomps = endPositionOfBraPair - startPositionOfBraPair;

        // compute arguments of Boys function values of maximum order

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to Boys function arguments

            int32_t ioff = posIntegralInBuffer + maxOrderOfIntegral * ncomps + i - startPositionOfBraPair;

            double* bvals = (double*)((char*)primBufferData + ioff * pitchOfBufferData);

            // set up pointers to Obara-Saika prefactors

            int32_t joff = 4 * (i - startPositionOfBraPair);

            const double* fz = (double*)((char*)osFactorsData + (joff  + 1) * pitchOfFactorsData);

            // set up pointers to R(PQ) distances

            int32_t koff = 3 * (i - startPositionOfBraPair);

            double* pqx = (double*)((char*)pqDistancesData + koff * pitchOfDistancesPQData);

            double* pqy = (double*)((char*)pqDistancesData + (koff + 1) * pitchOfDistancesPQData);

            double* pqz = (double*)((char*)pqDistancesData + (koff + 2) * pitchOfDistancesPQData);

            // compute Boys function of maximum order

            barg = fz[tid] * (pqx[tid] * pqx[tid] + pqy[tid] * pqy[tid] + pqz[tid] * pqz[tid]);

            bvals[tid] = boys (maxOrderOfIntegral, barg);
        }
    }
}

void
launchKernelEriForSSSS(      double* primBufferData,
                       const size_t  pitchOfBufferData,
                       const int32_t posIntegralInBuffer,
                       const int32_t maxOrderOfIntegral,
                       const double* osFactorsData,
                       const size_t  pitchOfFactorsData,
                       const double* pqDistancesData,
                       const size_t  pitchOfDistancesPQData,
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
    gpu::kernelEriForSSSS<<<gridSize, blockSize>>>(primBufferData, pitchOfBufferData,
                                                   posIntegralInBuffer, maxOrderOfIntegral,
                                                   osFactorsData, pitchOfFactorsData,
                                                   braGtoPairsData, pitchOfBraGtoPairsData,
                                                   ketGtoPairsData, pitchOfKetGtoPairsData,
                                                   startPositionOfBraPair, endPositionOfBraPair,
                                                   nKetPrimPairs);
}

}  // namespace gpu

