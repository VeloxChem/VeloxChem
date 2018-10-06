//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef KetHrrFunc_hpp
#define KetHrrFunc_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"

namespace kethrrfunc { // kethrrfunc namespace
    
    /**
     Computes batch of contracted (SX|g(r,r')|PP) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPP(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PD) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPD(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPF(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPG(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPH(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PI) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPI(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|PK) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXPK(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|DD) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXDD(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|DF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXDF(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|DG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXDG(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|DH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXDH(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|DI) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXDI(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|FF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXFF(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|FG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXFG(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|FH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXFH(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (SX|g(r,r')|GG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer for ket side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForSXGG(      CMemBlock2D<double>&  ketBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  cdDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetContrPairs,
                                      const int32_t               iContrPair);
    
} // kethrrfunc namespace

#endif /* KetHrrFunc_hpp */
