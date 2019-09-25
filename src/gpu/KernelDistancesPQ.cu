//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KernelDistancesPQ.hpp"

namespace gpu {  // gpu namespace

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
    // FIX ME: add CUDA kernel body 
}

}  // namespace gpu
