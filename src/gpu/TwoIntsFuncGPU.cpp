//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoIntsFuncGPU.hpp"

#include "CudaGenFunc.hpp"
#include "EriRecFuncGPU.hpp"

#ifdef ENABLE_GPU
#include "KernelDistancesPQ.hpp"
#include "KernelElectronRepulsionFactors.hpp"
#include "KernelCoordinatesW.hpp"
#include "KernelDistancesWP.hpp"
#include "KernelDistancesWQ.hpp"
#endif

namespace twointsgpu { // twointsgpu namespace
    
    void
    compDistancesPQ(      double*         pqDistancesData,
                    const size_t          pitchOfDistancesData,
                    const double*         braGtoPairsData,
                    const size_t          pitchOfBraGtoPairsData,
                    const double*         ketGtoPairsData,
                    const size_t          pitchOfKetGtoPairsData,
                    const CGtoPairsBlock& braGtoPairsBlock,
                    const int32_t         nKetPrimPairs,
                    const int32_t         iContrPair,
                    const CCudaDevices*   cudaDevices)
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
    compFactorsForElectronRepulsion(      double*         osFactorsData,
                                    const size_t          pitchOfFactorsData,
                                    const double*         braGtoPairsData,
                                    const size_t          pitchOfBraGtoPairsData,
                                    const double*         ketGtoPairsData,
                                    const size_t          pitchOfKetGtoPairsData,
                                    const CGtoPairsBlock& braGtoPairsBlock,
                                    const CGtoPairsBlock& ketGtoPairsBlock,
                                    const int32_t         nKetPrimPairs,
                                    const int32_t         iContrPair,
                                    const CCudaDevices*   cudaDevices)
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
    compCoordinatesW(      double*         wCoordinatesData,
                     const size_t          pitchOfCoordinatesData,
                     const double*         osFactorsData,
                     const size_t          pitchOfFactorsData,
                     const double*         braGtoPairsData,
                     const size_t          pitchOfBraGtoPairsData,
                     const double*         ketGtoPairsData,
                     const size_t          pitchOfKetGtoPairsData,
                     const CGtoPairsBlock& braGtoPairsBlock,
                     const int32_t         nKetPrimPairs,
                     const int32_t         iContrPair,
                     const CCudaDevices*   cudaDevices)
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
    
    void
    compDistancesWP(      double*         wpDistancesData,
                    const size_t          pitchOfDistancesData,
                    const double*         wCoordinatesData,
                    const size_t          pitchOfCoordinatesData,
                    const double*         braGtoPairsData,
                    const size_t          pitchOfBraGtoPairsData,
                    const CGtoPairsBlock& braGtoPairsBlock,
                    const int32_t         nKetPrimPairs,
                    const int32_t         iContrPair,
                    const CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
    // skip computation for zero angular momentum on bra side
        
    if ((braGtoPairsBlock.getBraAngularMomentum() == 0) &&
        (braGtoPairsBlock.getKetAngularMomentum() == 0)) return;
        
    // set up GTOs pair position on bra side
    
    auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
    
    auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
    
    //  determine execution grid on GPU device
    
    auto bsize = cudaDevices->getGridBlockSize();
    
    auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
    
    // execute CUDA kernel: R(WP) distances
    
    gpu::launchKernelForDistancesWP(wpDistancesData, pitchOfDistancesData, wCoordinatesData, pitchOfCoordinatesData,
                                    braGtoPairsData, pitchOfBraGtoPairsData, spos, epos, nKetPrimPairs, gsize, bsize);
#endif
    }
    
    void
    compDistancesWQ(      double*         wqDistancesData,
                    const size_t          pitchOfDistancesData,
                    const double*         wCoordinatesData,
                    const size_t          pitchOfCoordinatesData,
                    const double*         ketGtoPairsData,
                    const size_t          pitchOfKetGtoPairsData,
                    const CGtoPairsBlock& braGtoPairsBlock,
                    const CGtoPairsBlock& ketGtoPairsBlock,
                    const int32_t         nKetPrimPairs,
                    const int32_t         iContrPair,
                    const CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
        // skip computation for zero angular momentum on bra side
        
        if ((ketGtoPairsBlock.getBraAngularMomentum() == 0) &&
            (ketGtoPairsBlock.getKetAngularMomentum() == 0)) return;
        
        // set up GTOs pair position on bra side
        
        auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
        
        auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
        
        //  determine execution grid on GPU device
        
        auto bsize = cudaDevices->getGridBlockSize();
        
        auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
        
        // execute CUDA kernel: R(WQ) distances
        
        gpu::launchKernelForDistancesWQ(wqDistancesData, pitchOfDistancesData, wCoordinatesData, pitchOfCoordinatesData,
                                        ketGtoPairsData, pitchOfKetGtoPairsData, spos, epos, nKetPrimPairs, gsize, bsize);
#endif
    }
    
    void
    compPrimElectronRepulsionIntsOnGPU(      double*         primBufferData,
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
                                       const CCudaDevices*   cudaDevices)
    {
#ifdef ENABLE_GPU
        // set up angular momentum
        
        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();
        
        auto kang = ketGtoPairsBlock.getBraAngularMomentum() + ketGtoPairsBlock.getKetAngularMomentum();
        
        // set up GTOs pair position on bra side
        
        auto spos = (braGtoPairsBlock.getStartPositions())[iContrPair];
        
        auto epos = (braGtoPairsBlock.getEndPositions())[iContrPair];
        
        //  determine execution grid on GPU device
        
        auto bsize = cudaDevices->getGridBlockSize();
        
        auto gsize = gpu::getNumberOfGridBlocks(nKetPrimPairs, bsize);
        
        // primitive integrals (SS||SS) integrals
        
        if ((bang == 0) && (kang == 0))
        {
            erirecgpu::compElectronRepulsionForSSSS(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SS||SP) integrals
        
        if ((bang == 0) && (kang == 1))
        {
            erirecgpu::compElectronRepulsionForSSSP(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wqDistancesData, pitchOfDistancesWQData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SP||SS) integrals
        
        if ((bang == 1) && (kang == 0))
        {
            erirecgpu::compElectronRepulsionForSPSS(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wpDistancesData, pitchOfDistancesWPData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SS||SD) integrals
        
        if ((bang == 0) && (kang == 2))
        {
            erirecgpu::compElectronRepulsionForSSSD(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wqDistancesData, pitchOfDistancesWQData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SD||SS) integrals
        
        if ((bang == 2) && (kang == 0))
        {
            erirecgpu::compElectronRepulsionForSDSS(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wpDistancesData, pitchOfDistancesWPData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SS||SF) integrals
        
        if ((bang == 0) && (kang == 3))
        {
            erirecgpu::compElectronRepulsionForSSSF(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wqDistancesData, pitchOfDistancesWQData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SF||SS) integrals
        
        if ((bang == 3) && (kang == 0))
        {
            erirecgpu::compElectronRepulsionForSFSS(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wpDistancesData, pitchOfDistancesWPData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SS||SG) integrals
        
        if ((bang == 0) && (kang == 4))
        {
            erirecgpu::compElectronRepulsionForSSSG(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wqDistancesData, pitchOfDistancesWQData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
        
        // primitive integrals (SG||SS) integrals
        
        if ((bang == 4) && (kang == 0))
        {
            erirecgpu::compElectronRepulsionForSGSS(primBufferData, pitchOfBufferData,
                                                    osFactorsData, pitchOfFactorsData,
                                                    pqDistancesData, pitchOfDistancesPQData,
                                                    wpDistancesData, pitchOfDistancesWPData,
                                                    braGtoPairsData, pitchOfBraGtoPairsData,
                                                    ketGtoPairsData, pitchOfKetGtoPairsData,
                                                    spos, epos, nKetPrimPairs, gsize, bsize);
            
            return;
        }
#endif
    }
    
} // intsfunc namespace
