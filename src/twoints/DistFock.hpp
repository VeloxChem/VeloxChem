//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DistFock_hpp
#define DistFock_hpp

#include <cstdlib>
#include <cstdint>

#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"
#include "FockContainer.hpp"

namespace distfock { // distfock namespace
    
    /**
     Distributes vector of spherical integrals into spin restricted closed shell
     Hatree-Fock matrix (2J + K)

     @param fockContainer the Fock container.
     @param iFockMatrix the index of Fock submatrix in Fock container.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param spherInts the spherical integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distRestJK(      CFockContainer&      fockContainer,
                    const int32_t              iFockMatrix,
                    const double*              densityMatrix,
                    const int32_t              nDensityColumns,
                    const CMemBlock2D<double>& spherInts,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              nKetContrPairs,
                    const int32_t              iContrPair);
    
    /**
     Distributes vector of spherical integrals into spin restricted closed shell
     Hatree-Fock matrix (J)
     
     @param fockContainer the Fock container.
     @param iFockMatrix the index of Fock submatrix in Fock container.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param spherInts the spherical integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distRestJ(      CFockContainer&      fockContainer,
                   const int32_t              iFockMatrix,
                   const double*              densityMatrix,
                   const int32_t              nDensityColumns,
                   const CMemBlock2D<double>& spherInts,
                   const CGtoPairsBlock&      braGtoPairsBlock,
                   const CGtoPairsBlock&      ketGtoPairsBlock,
                   const int32_t              nKetContrPairs,
                   const int32_t              iContrPair);
    
    /**
     Distributes vector of spherical integrals into spin restricted closed shell
     Hatree-Fock matrix (K)
     
     @param fockContainer the Fock container.
     @param iFockMatrix the index of Fock submatrix in Fock container.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param spherInts the spherical integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distRestK(      CFockContainer&      fockContainer,
                   const int32_t              iFockMatrix,
                   const double*              densityMatrix,
                   const int32_t              nDensityColumns,
                   const CMemBlock2D<double>& spherInts,
                   const CGtoPairsBlock&      braGtoPairsBlock,
                   const CGtoPairsBlock&      ketGtoPairsBlock,
                   const int32_t              nKetContrPairs,
                   const int32_t              iContrPair);
    
    /**
     Distributes vector of spherical integrals into spin restricted closed shell
     general Coulomb matrix (J)
     
     @param fockContainer the Fock container.
     @param iFockMatrix the index of Fock submatrix in Fock container.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param spherInts the spherical integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distRestGenJ(      CFockContainer&      fockContainer,
                      const int32_t              iFockMatrix,
                      const double*              densityMatrix,
                      const int32_t              nDensityColumns,
                      const CMemBlock2D<double>& spherInts,
                      const CGtoPairsBlock&      braGtoPairsBlock,
                      const CGtoPairsBlock&      ketGtoPairsBlock,
                      const int32_t              nKetContrPairs,
                      const int32_t              iContrPair);
    
    /**
     Distributes vector of spherical integrals into spin restricted closed shell
     general exchange matrix (K)
     
     @param fockContainer the Fock container.
     @param iFockMatrix the index of Fock submatrix in Fock container.
     @param densityMatrix the constant pointer to AO density matrix.
     @param nDensityColumns the number of columns in AO density matrix.
     @param spherInts the spherical integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distRestGenK(      CFockContainer&      fockContainer,
                      const int32_t              iFockMatrix,
                      const double*              densityMatrix,
                      const int32_t              nDensityColumns,
                      const CMemBlock2D<double>& spherInts,
                      const CGtoPairsBlock&      braGtoPairsBlock,
                      const CGtoPairsBlock&      ketGtoPairsBlock,
                      const int32_t              nKetContrPairs,
                      const int32_t              iContrPair);
    
    /**
    Distributes vector of spherical integrals into spin restricted closed shell
    general Fock matrix (2J - K)
    
    @param fockContainer the Fock container.
    @param iFockMatrix the index of Fock submatrix in Fock container.
    @param densityMatrix the constant pointer to AO density matrix.
    @param nDensityColumns the number of columns in AO density matrix.
    @param spherInts the spherical integrals buffer.
    @param braGtoPairsBlock the GTOs pairs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contracted GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void distRestGenJK(      CFockContainer&      fockContainer,
                       const int32_t              iFockMatrix,
                       const double*              densityMatrix,
                       const int32_t              nDensityColumns,
                       const CMemBlock2D<double>& spherInts,
                       const CGtoPairsBlock&      braGtoPairsBlock,
                       const CGtoPairsBlock&      ketGtoPairsBlock,
                       const int32_t              nKetContrPairs,
                       const int32_t              iContrPair);
    
} // distfock namespace

#endif /* DistFock_hpp */
