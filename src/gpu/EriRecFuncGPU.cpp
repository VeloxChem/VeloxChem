//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "EriRecFuncGPU.hpp"

#ifdef ENABLE_GPU
#include "KernelEriForSSSS.hpp"
#endif

namespace erirecgpu { // erirecgpu namespace
   
    void
    compElectronRepulsionForSSSS(      double* primBufferData,
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
#ifdef ENABLE_GPU
        gpu::launchKernelEriForSSSS(primBufferData, pitchOfBufferData, posIntegralInBuffer, maxOrderOfIntegral,
                                    osFactorsData, pitchOfFactorsData, pqDistancesData, pitchOfDistancesPQData,
                                    braGtoPairsData, pitchOfBraGtoPairsData, ketGtoPairsData, pitchOfKetGtoPairsData,
                                    startPositionOfBraPair, endPositionOfBraPair, nKetPrimPairs, gridSize, blockSize);
#endif 
    }
    
} // erirecgpu namespace
