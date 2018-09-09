//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef BraHrrFunc_hpp
#define BraHrrFunc_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"

namespace brahrrfunc { // brahrrfunc namespace
    
    /**
     Computes batch of contracted (PP|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPPXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PD|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPDXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PF|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPFXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PG|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPGXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PH|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPHXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PI|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPIXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (PK|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForPKXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (DD|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForDDXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (DF|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForDFXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (DG|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForDGXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (DH|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForDHXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (DI|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForDIXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (FF|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForFFXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (FG|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForFGXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (FH|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForFHXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of contracted (GG|g(r,r')|XX) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param braBuffer the horizontal recursion buffer for bra side.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compElectronRepulsionForGGXX(      CMemBlock2D<double>&  braBuffer,
                                      const CVecFourIndexes&      recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
} // brahrrfunc namespace

#endif /* BraHrrFunc_hpp */
