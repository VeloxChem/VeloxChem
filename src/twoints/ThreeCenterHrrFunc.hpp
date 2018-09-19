//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ThreeCenterHrrFunc_hpp
#define ThreeCenterHrrFunc_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"

namespace t3hrrfunc { // t3hrrfunc namespace
 
    /**
     Computes batch of half transformed (X|g(r,r')|PP) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPP(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PD) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPD(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPF(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPG(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPH(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PI) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPI(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|PK) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPK(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|DD) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXDD(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|DF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXDF(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|DG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXDG(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|DH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXDH(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|DI) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXDI(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|FF) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXFF(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|FG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXFG(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|FH) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXFH(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    /**
     Computes batch of half transformed (X|g(r,r')|GG) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXGG(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
    
} // t3hrrfunc namespace

#endif /* ThreeCenterHrrFunc_hpp */
