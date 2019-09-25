//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DeviceProp_hpp
#define DeviceProp_hpp

#include <cstdint>

namespace gpu {

    /**
     Launches CUDA kernel for computation of R(PQ) = P - Q distances.

     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(P-Q) data on CUDA compute capable device.
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
    void launchKernelForDistancesPQ(      double* pqDistancesData,
                                    const size_t  pitchOfDistancesData,
                                    const double* braGtoPairsData,
                                    const size_t  pitchOfBraGtoPairsData,
                                    const double* ketGtoPairsData,
                                    const size_t  pitchOfKetGtoPairsData,
                                    const int32_t startPositionOfBraPair,
                                    const int32_t endPositionOfBraPair,
                                    const int32_t nKetPrimPairs,
                                    const int32_t gridSize,
                                    const int32_t blockSize);
    
    /**
     CUDA kernel for computation of R(PQ) = P - Q distances.

     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(P-Q) data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param startPositionOfBraPair the start position of GTOs pair on bra side.
     @param endPositionOfBraPair the end position of GTOs pair on bra side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     */
    __global__ void kernelDistancesPQ(      double* pqDistancesData,
                                      const size_t  pitchOfDistancesData,
                                      const double* braGtoPairsData,
                                      const size_t  pitchOfBraGtoPairsData,
                                      const double* ketGtoPairsData,
                                      const size_t  pitchOfKetGtoPairsData,
                                      const int32_t startPositionOfBraPair,
                                      const int32_t endPositionOfBraPair,
                                      const int32_t nKetPrimPairs);
}

#endif
