//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef KernelDistancesWP_hpp
#define KernelDistancesWP_hpp

#include <cstdint>

namespace gpu {
    /**
     Launches CUDA kernel for computation of combined distances R(WP) = W - P.

     @param wpDistancesData the vector of distances R(WP) = W - P on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of distances R(WP) = W - P on CUDA compute capable device.
     @param wCoordinatesData the vector of combined coordinates W on CUDA compute capable device.
     @param pitchOfCoordinatesData the pitch of combined coordinates W on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param startPositionOfBraPair the start position of GTOs pair on bra side.
     @param endPositionOfBraPair the end position of GTOs pair on bra side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param gridSize the execution grid size on device.
     @param blockSize the execution block size on device.
     */
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
                                    const int32_t blockSize);
}

#endif
