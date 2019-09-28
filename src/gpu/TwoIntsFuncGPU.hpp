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

     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
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
    void compFactorsForElectronRepulsion(        double*         osFactorsData,
                                         const   size_t          pitchOfFactorsData,
                                         const   double*         braGtoPairsData,
                                         const   size_t          pitchOfBraGtoPairsData,
                                         const   double*         ketGtoPairsData,
                                         const   size_t          pitchOfKetGtoPairsData,
                                         const   CGtoPairsBlock& braGtoPairsBlock,
                                         const   CGtoPairsBlock& ketGtoPairsBlock,
                                         const   int32_t         nKetPrimPairs,
                                         const   int32_t         iContrPair,
                                         const   CCudaDevices*   cudaDevices);
    
    /**
     Computes coordinates of combined bra and ket GTOs pairs.
     
     @param wCoordinatesData the vector of combined W coordinates on CUDA compute capable device.
     @param pitchOfCoordinatesData the pitch of W coordinates data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param braGtoPairsBlock the bra GTOs pairs block object.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param cudaDevices the CUDA compute capable devices.
     */
    void compCoordinatesW(        double*         wCoordinatesData,
                          const   size_t          pitchOfCoordinatesData,
                          const   double*         osFactorsData,
                          const   size_t          pitchOfFactorsData,
                          const   double*         braGtoPairsData,
                          const   size_t          pitchOfBraGtoPairsData,
                          const   double*         ketGtoPairsData,
                          const   size_t          pitchOfKetGtoPairsData,
                          const   CGtoPairsBlock& braGtoPairsBlock,
                          const   int32_t         nKetPrimPairs,
                          const   int32_t         iContrPair,
                          const   CCudaDevices*   cudaDevices);
    
    /**
     Computes distances R(WP) = W - P.
     
     @param wpDistancesData the vector of R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(WP) distances data on CUDA compute capable device.
     @param wCoordinatesData the vector of combined W coordinates on CUDA compute capable device.
     @param pitchOfCoordinatesData the pitch of W coordinates data on CUDA compute capable device.
     @param braGtoPairsData the bra GTOs pairs data on on CUDA compute capable device.
     @param pitchOfBraGtoPairsData the pitch of bra GTOs pairs data on on CUDA compute capable device.
     @param braGtoPairsBlock the bra GTOs pairs block object.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param cudaDevices the CUDA compute capable devices.
     */
    void compDistancesWP(        double*         wpDistancesData,
                         const   size_t          pitchOfDistancesData,
                         const   double*         wCoordinatesData,
                         const   size_t          pitchOfCoordinatesData,
                         const   double*         braGtoPairsData,
                         const   size_t          pitchOfBraGtoPairsData,
                         const   CGtoPairsBlock& braGtoPairsBlock,
                         const   int32_t         nKetPrimPairs,
                         const   int32_t         iContrPair,
                         const   CCudaDevices*   cudaDevices);
    
    /**
     Computes distances R(WQ) = W - Q.
     
     @param wqDistancesData the vector of R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesData the pitch of R(WQ) distances data on CUDA compute capable device.
     @param wCoordinatesData the vector of combined W coordinates on CUDA compute capable device.
     @param pitchOfCoordinatesData the pitch of W coordinates data on CUDA compute capable device.
     @param ketGtoPairsData the ket GTOs pairs data on on CUDA compute capable device.
     @param pitchOfKetGtoPairsData the pitch of ket GTOs pairs data on on CUDA compute capable device.
     @param braGtoPairsBlock the bra GTOs pairs block object.
     @param ketGtoPairsBlock the ket GTOs pairs block object.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param cudaDevices the CUDA compute capable devices.
     */
    void compDistancesWQ(        double*         wqDistancesData,
                         const   size_t          pitchOfDistancesData,
                         const   double*         wCoordinatesData,
                         const   size_t          pitchOfCoordinatesData,
                         const   double*         ketGtoPairsData,
                         const   size_t          pitchOfKetGtoPairsData,
                         const   CGtoPairsBlock& braGtoPairsBlock,
                         const   CGtoPairsBlock& ketGtoPairsBlock,
                         const   int32_t         nKetPrimPairs,
                         const   int32_t         iContrPair,
                         const   CCudaDevices*   cudaDevices);
    
    /**
     Computes primitive electron repulsion integrals using vertical Obara-Saika recursion.

     @param primBufferData the vertical recursion buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of vertical recursion buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData vector of R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(WP) distances data on CUDA compute capable device.
     @param wqDistancesData the vector of R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(WQ) distances data on CUDA compute capable device.
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
    void compPrimElectronRepulsionIntsOnGPU(      double*         primBufferData,
                                            const size_t          pitchOfBufferData,
                                            const double*         osFactorsData,
                                            const size_t          pitchOfFactorsData,
                                            const double*         pqDistancesData,
                                            const size_t          pitchOfDistancesPQData,
                                            const double*         wpDistancesData,
                                            const size_t          pitchOfDistancesWPData,
                                            const double*         wqDistancesData,
                                            const size_t          pitchOfDistancesWQData,
                                            const double*         braGtoPairsData,
                                            const size_t          pitchOfBraGtoPairsData,
                                            const double*         ketGtoPairsData,
                                            const size_t          pitchOfKetGtoPairsData,
                                            const CGtoPairsBlock& braGtoPairsBlock,
                                            const CGtoPairsBlock& ketGtoPairsBlock,
                                            const int32_t         nKetPrimPairs,
                                            const int32_t         iContrPair,
                                            const CCudaDevices*   cudaDevices);
    
} // twointsgpu namespace

#endif /* TwoIntsFuncGPU_hpp */
