//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TwoIntsFuncGPU_hpp
#define TwoIntsFuncGPU_hpp

#include <cstdint>

#include "GtoPairsBlock.hpp"
#include "CudaDevices.hpp"

namespace twointsgpu { // twointsgpu namespace

    /**
     Computes distances between combined centers of contracted GTOs pair on bra
     side with all combined centers of contracted GTOs pairs on ket side.

     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(P-Q) data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param braGtoPairsBlock the bra GTOs pairs block object.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param cudaDevices the CUDA compute capable devices.
     */
    void compDistancesPQ(        double*         pqDistancesData,
                         const   size_t          pitchOfDistancesData,
                         const   double*         braGtoPairsData,
                         const   size_t          pitchOfBraGtoPairsData,
                         const   double*         ketGtoPairsData,
                         const   size_t          pitchOfKetGtoPairsData,
                         const   CGtoPairsBlock& braGtoPairsBlock,
                         const   int32_t         nKetPrimPairs,
                         const   int32_t         iContrPair,
                         const   CCudaDevices*   cudaDevices);
    
    /**
     Computes primitive Obara-Saika recursion factors for electron repulsion integrals.

     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(P-Q) data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param braGtoPairsBlock the bra GTOs pairs block object.
     @param ketGtoPairsBlock the ket GTOs pairs block object.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param cudaDevices the CUDA compute capable devices.
     */
    void compFactorsForElectronRepulsion(        double*         pqDistancesData,
                                         const   size_t          pitchOfDistancesData,
                                         const   double*         braGtoPairsData,
                                         const   size_t          pitchOfBraGtoPairsData,
                                         const   double*         ketGtoPairsData,
                                         const   size_t          pitchOfKetGtoPairsData,
                                         const   CGtoPairsBlock& braGtoPairsBlock,
                                         const   CGtoPairsBlock& ketGtoPairsBlock,
                                         const   int32_t         nKetPrimPairs,
                                         const   int32_t         iContrPair,
                                         const   CCudaDevices*   cudaDevices);
    
} // twointsgpu namespace

#endif /* TwoIntsFuncGPU_hpp */
