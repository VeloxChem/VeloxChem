//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef EriRecFuncGPU_hpp
#define EriRecFuncGPU_hpp

#include <cstdint>

#include "GtoPairsBlock.hpp"
#include "CudaDevices.hpp"

namespace erirecgpu { // erirecgpu namespace

    /**
     Computes primitive (SS||SS)^(m) integrals.

     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param posIntegralInBuffer the position of (SS||SS)^(m) in primitives buffer.
     @param maxOrderOfIntegral the maximum order of (SS||SS)^(m) integrals.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param startPositionOfBraPair the start position of GTOs pair on bra side.
     @param endPositionOfBraPair the end position of GTOs pair on bra side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param gridSize the execution grid size on device.
     @param blockSize the execution block size on device.
     */
    void compElectronRepulsionForSSSS(      double* primBufferData,
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
                                      const int32_t blockSize);
    
} // erirecgpu namespace

#endif /* EriRecFuncGPU_hpp */
