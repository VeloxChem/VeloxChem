//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoIntsFuncGPU.hpp"

#include "CudaGenFunc.hpp"

#ifdef ENABLE_GPU
#include "KernelDistancesPQ.hpp"
#include "KernelElectronRepulsionFactors.hpp"
#include "KernelCoordinatesW.hpp"
#endif

namespace twointsgpu { // twointsgpu namespace
    
    void
    compDistancesPQ(        double*         pqDistancesData,
                    const   size_t          pitchOfDistancesData,
                    const   double*         braGtoPairsData,
                    const   size_t          pitchOfBraGtoPairsData,
                    const   double*         ketGtoPairsData,
                    const   size_t          pitchOfKetGtoPairsData,
                    const   CGtoPairsBlock& braGtoPairsBlock,
                    const   int32_t         nKetPrimPairs,
                    const   int32_t         iContrPair,
                    const   CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
        // set up GTOs pair position on bra side
        
        auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
        
        auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
        
        //  determine execution grid on GPU device
        
        auto bsize = cudaDevices->getGridBlockSize();
        
        auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
        
        // execute CUDA kernel: R(PQ) = P - Q
        
        gpu::launchKernelForDistancesPQ(pqDistancesData, pitchOfDistancesData, braGtoPairsData, pitchOfBraGtoPairsData,
                                        ketGtoPairsData, pitchOfKetGtoPairsData, spos, epos, nKetPrimPairs,
                                        gsize, bsize); 
#endif
    }
    
    void
    compFactorsForElectronRepulsion(        double*         osFactorsData,
                                    const   size_t          pitchOfFactorsData,
                                    const   double*         braGtoPairsData,
                                    const   size_t          pitchOfBraGtoPairsData,
                                    const   double*         ketGtoPairsData,
                                    const   size_t          pitchOfKetGtoPairsData,
                                    const   CGtoPairsBlock& braGtoPairsBlock,
                                    const   CGtoPairsBlock& ketGtoPairsBlock,
                                    const   int32_t         nKetPrimPairs,
                                    const   int32_t         iContrPair,
                                    const   CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
        // set up angular momentum data
        
        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();
        
        auto kang = ketGtoPairsBlock.getBraAngularMomentum() + ketGtoPairsBlock.getKetAngularMomentum();
        
        // set up GTOs pair position on bra side
        
        auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
        
        auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
        
        //  determine execution grid on GPU device
        
        auto bsize = cudaDevices->getGridBlockSize();
        
        auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
        
        // execute CUDA kernel: Obara-Saika recursion factors for electron repulsion integrals
        
        gpu::launchKernelForElectronRepulsionFactors(osFactorsData, pitchOfFactorsData, braGtoPairsData,
                                                     pitchOfBraGtoPairsData, ketGtoPairsData, pitchOfKetGtoPairsData,
                                                     bang, kang, spos, epos, nKetPrimPairs, gsize, bsize);
#endif
    }
    
    void
    compCoordinatesW(        double*         wCoordinatesData,
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
                     const   CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
        // set up GTOs pair position on bra side
        
        auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
        
        auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
        
        //  determine execution grid on GPU device
        
        auto bsize = cudaDevices->getGridBlockSize();
        
        auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
        
        // execute CUDA kernel: W coordinates
        
        gpu::launchKernelForCoordinatesW(wCoordinatesData, pitchOfCoordinatesData, osFactorsData, pitchOfFactorsData,
                                         braGtoPairsData, pitchOfBraGtoPairsData, ketGtoPairsData, pitchOfKetGtoPairsData,
                                         spos, epos, nKetPrimPairs, gsize, bsize);
#endif
    }
    
} // intsfunc namespace
