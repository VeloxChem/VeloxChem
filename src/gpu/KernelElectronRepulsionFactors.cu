//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelDistancesPQ.hpp"

namespace gpu {  // gpu namespace

__global__ void kernelElectronRepulsionFactors(      double* osFactorsData,
                                               const size_t  pitchOfFactorsData,
                                               const double* braGtoPairsData,
                                               const size_t  pitchOfBraGtoPairsData,
                                               const double* ketGtoPairsData,
                                               const size_t  pitchOfKetGtoPairsData,
                                               const int32_t braAngularMomentum,
                                               const int32_t ketAngularMomentum,
                                               const int32_t startPositionOfBraPair,
                                               const int32_t endPositionOfBraPair,
                                               const int32_t nKetPrimPairs)
{
    int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < nKetPrimPairs)
    {
        // set up pointers to primitve pair factors on bra side

        double* bfxi = braGtoPairsData;

        double* boxi = (double*)((char*)braGtoPairsData + pitchOfBraGtoPairsData);

        // set up pointers to primitve pair factors on ket side

        double* kfxi = ketGtoPairsData;

        double* koxi = (double*)((char*)ketGtoPairsData + pitchOfKetGtoPairsData);

        // compute Obara-Saika factors

        for (int32_t i = startPositionOfBraPair; i < endPositionOfBraPair; i++)
        {
            // set up pointers to Obara-Saika factors 

            int32_t ioff = 4 * (i - startPositionOfBraPair);

            // auto fb = bfxi[i];

            double* fx = (double*)((char*)osFactorsData + ioff * pitchOfFactorsData);

            double* fz = (double*)((char*)osFactorsData + (ioff + 1) * pitchOfFactorsData);

            // compute fx and fz factors

            fx[tid] = 1.0 / (bfxi[i] + kfxi[tid]);

            fz[tid] = bfxi[i] * kfxi[tid] * fx[tid];

            if (braAngularMomentum > 1)
            {
                double* ta = (double*)((char*)osFactorsData + (ioff + 2) * pitchOfFactorsData);

                ta[tid] = boxi[i] * fz[tid];
            }

            if (ketAngularMomentum > 1)
            {
                double* td = (double*)((char*)osFactorsData + (ioff + 3) * pitchOfFactorsData);

                td[tid] = koxi[tid] * fz[tid];
            }
        }
    }
}

void
launchKernelForElectronRepulsionFactors(      double* osFactorsData,
                                        const size_t  pitchOfFactorsData,
                                        const double* braGtoPairsData,
                                        const size_t  pitchOfBraGtoPairsData,
                                        const double* ketGtoPairsData,
                                        const size_t  pitchOfKetGtoPairsData,
                                        const int32_t braAngularMomentum,
                                        const int32_t ketAngularMomentum,
                                        const int32_t startPositionOfBraPair,
                                        const int32_t endPositionOfBraPair,
                                        const int32_t nKetPrimPairs,
                                        const int32_t gridSize,
                                        const int32_t blockSize)
{
    gpu::kernelElectronRepulsionFactors<<<gridSize, blockSize>>>(osFactorsData, pitchOfFactorsData,
                                                                 braGtoPairsData, pitchOfBraGtoPairsData,
                                                                 ketGtoPairsData, pitchOfKetGtoPairsData,
                                                                 braAngularMomentum, ketAngularMomentum,
                                                                 startPositionOfBraPair, endPositionOfBraPair,
                                                                 nKetPrimPairs);
}

}  // namespace gpu

