//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DistMaxDensity_hpp
#define DistMaxDensity_hpp

#include <cstdint>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"

namespace distmaxden { // distmaxden namespace
    
    /**
     Determines maximum density elements for vector of GTOs pairs on ket side
     and GTOs pair on bra side in case of spin restricted Fock matrix, 2J + K.

     @param maxDensityElements maxDensityElements the vector of maximum density
            elements.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxRestDenJK(      CMemBlock<double>&   maxDensityElements,
                         const double*              densityMatrix,
                         const int32_t              nDensityColumns,
                         const CGtoPairsBlock&      braGtoPairsBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              nKetContrPairs,
                         const int32_t              iContrPair);
    
    /**
     Determines maximum density elements for vector of GTOs pairs on ket side
     and GTOs pair on bra side in case of spin restricted Fock matrix, J.
     
     @param maxDensityElements maxDensityElements the vector of maximum density
            elements.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxRestDenJ(      CMemBlock<double>&   maxDensityElements,
                        const double*              densityMatrix,
                        const int32_t              nDensityColumns,
                        const CGtoPairsBlock&      braGtoPairsBlock,
                        const CGtoPairsBlock&      ketGtoPairsBlock,
                        const int32_t              nKetContrPairs,
                        const int32_t              iContrPair);
    
    /**
     Determines maximum density elements for vector of GTOs pairs on ket side
     and GTOs pair on bra side in case of spin restricted Fock matrix, K.
     
     @param maxDensityElements maxDensityElements the vector of maximum density
     elements.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxRestDenK(      CMemBlock<double>&   maxDensityElements,
                        const double*              densityMatrix,
                        const int32_t              nDensityColumns,
                        const CGtoPairsBlock&      braGtoPairsBlock,
                        const CGtoPairsBlock&      ketGtoPairsBlock,
                        const int32_t              nKetContrPairs,
                        const int32_t              iContrPair);
    
    /**
     Determines maximum density elements for vector of GTOs pairs on ket side
     and GTOs pair on bra side in case of spin restricted general Coulomb
     matrix, J.
     
     @param maxDensityElements maxDensityElements the vector of maximum density
     elements.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxRestGenDenJ(      CMemBlock<double>&   maxDensityElements,
                           const double*              densityMatrix,
                           const int32_t              nDensityColumns,
                           const CGtoPairsBlock&      braGtoPairsBlock,
                           const CGtoPairsBlock&      ketGtoPairsBlock,
                           const int32_t              nKetContrPairs,
                           const int32_t              iContrPair);
    
    /**
     Determines maximum density elements for vector of GTOs pairs on ket side
     and GTOs pair on bra side in case of spin restricted general exchange
     matrix, K.
     
     @param maxDensityElements maxDensityElements the vector of maximum density
     elements.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxRestGenDenK(      CMemBlock<double>&   maxDensityElements,
                           const double*              densityMatrix,
                           const int32_t              nDensityColumns,
                           const CGtoPairsBlock&      braGtoPairsBlock,
                           const CGtoPairsBlock&      ketGtoPairsBlock,
                           const int32_t              nKetContrPairs,
                           const int32_t              iContrPair);

} // distmaxden namespace

#endif /* DistMaxDensity_hpp */
