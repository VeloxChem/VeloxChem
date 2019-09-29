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
    
    /**
     Computes primitive (SS||SP)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSSSP(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SP||SS)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
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
    void compElectronRepulsionForSPSS(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
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
     Computes primitive (SS||SD)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSSSD(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SD||SS)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
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
    void compElectronRepulsionForSDSS(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
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
     Computes primitive (SS||SF)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSSSF(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SF||SS)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
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
    void compElectronRepulsionForSFSS(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
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
     Computes primitive (SS||SG)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSSSG(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SG||SS)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
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
    void compElectronRepulsionForSGSS(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
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
     Computes primitive (SP||SP)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSPSP(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SP||SD)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSPSD(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SD||SP)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSDSP(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SP||SF)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSPSF(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SF||SP)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSFSP(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SP||SG)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSPSG(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SG||SP)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSGSP(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SD||SD)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSDSD(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SD||SF)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSDSF(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SF||SD)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSFSD(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SD||SG)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSDSG(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SG||SD)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSGSD(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SF||SF)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSFSF(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SF||SG)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSFSG(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SG||SF)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSGSF(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
     Computes primitive (SG||SG)^(m) integrals.
     
     @param primBufferData the primitives buffer on CUDA compute capable device.
     @param pitchOfBufferData the pitch of primitives buffer data on CUDA compute capable device.
     @param osFactorsData the vector of Obara-Saika factors on CUDA compute capable device.
     @param pitchOfFactorsData the pitch of Obara-Saika factors data on CUDA compute capable device.
     @param pqDistancesData the vector of Cartesian R(PQ) = P - Q distances on CUDA compute capable device.
     @param pitchOfDistancesPQData the pitch of R(P-Q) data on CUDA compute capable device.
     @param wpDistancesData the vector of Cartesian R(WP) = W - P distances on CUDA compute capable device.
     @param pitchOfDistancesWPData the pitch of R(W-P) data on CUDA compute capable device.
     @param wqDistancesData the vector of Cartesian R(WQ) = W - Q distances on CUDA compute capable device.
     @param pitchOfDistancesWQData the pitch of R(W-Q) data on CUDA compute capable device.
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
    void compElectronRepulsionForSGSG(      double* primBufferData,
                                      const size_t  pitchOfBufferData,
                                      const double* osFactorsData,
                                      const size_t  pitchOfFactorsData,
                                      const double* pqDistancesData,
                                      const size_t  pitchOfDistancesPQData,
                                      const double* wpDistancesData,
                                      const size_t  pitchOfDistancesWPData,
                                      const double* wqDistancesData,
                                      const size_t  pitchOfDistancesWQData,
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
