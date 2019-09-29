//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef KernelDistancesWQ_hpp
#define KernelDistancesWQ_hpp

#include <cstdint>

namespace gpu {
    /**
     Launches CUDA kernel for computation of combined distances R(WQ) = W - Q.

     @param wqDistancesData the vector of distances R(WQ) = W - Q on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of distances R(WQ) = W - Q on CUDA compute capable device.
     @param wCoordinatesData the vector of combined coordinates W on CUDA compute capable device.
     @param pitchOfCoordinatesData the pitch of combined coordinates W on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param startPositionOfBraPair the start position of GTOs pair on bra side.
     @param endPositionOfBraPair the end position of GTOs pair on bra side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param gridSize the execution grid size on device.
     @param blockSize the execution block size on device.
     */
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
                                    const int32_t blockSize);
}

#endif
