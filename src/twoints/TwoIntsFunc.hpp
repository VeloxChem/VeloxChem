//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TwoIntsFunc_hpp
#define TwoIntsFunc_hpp

#include "MemBlock2D.hpp"
#include "MemBlock.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"

namespace twointsfunc { // twointsfunc namespace
    
    /**
     Computes distances between specific contracted GTO on bra side with all
     combined centers of contracted GTOs pairs on ket side.
     
     @param aqDistances the vector of Cartesian R(AQ) = A - Q distances.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesAQ(      CMemBlock2D<double>& aqDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for three center electron repulsion integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForThreeCenterElectronRepulsion(      CMemBlock2D<double>& osFactors,
                                                    const CGtoBlock&           braGtoBlock,
                                                    const CGtoPairsBlock&      ketGtoPairsBlock,
                                                    const int32_t              iContrGto);
    
    /**
     Computes coordinates of combined Gaussian functions, which is obtained by
     applying Gaussian product rule to thee Gaussian functions.
     
     @param wCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimension of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compCoordinatesForW(      CMemBlock2D<double>& wCoordinates,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nFactors,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoPairsBlock&      ketGtoPairsBlock,
                             const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center W of combined primitive GTOs
     and center A of primitive GTO on bra side.
     
     @param waDistances the vector of Cartesian R(WA) = W - A distances.
     @param wCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesWA(      CMemBlock2D<double>& waDistances,
                         const CMemBlock2D<double>& wCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center W of combined primitive GTOs
     and center D of primitive GTOs pair on ket side.
     
     @param wdDistances the vector of Cartesian R(WD) = W - D distances.
     @param wCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesWD(      CMemBlock2D<double>& wdDistances,
                         const CMemBlock2D<double>& wCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center W of combined primitive GTOs
     and combined center Q of primitive GTOs pair on ket side.
     
     @param wqDistances the vector of Cartesian R(WQ) = W - Q distances.
     @param wCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesWQ(      CMemBlock2D<double>& wqDistances,
                         const CMemBlock2D<double>& wCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes distances between combined centers of contracted GTOs pair on bra
     side with all combined centers of contracted GTOs pairs on ket side.
     
     @param pqDistances the vector of Cartesian R(PQ) = P - Q distances.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compDistancesPQ(      CMemBlock2D<double>& pqDistances,
                         const CGtoPairsBlock&      braGtoPairsBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              nKetPrimPairs,
                         const int32_t              iContrPair);
    
    /**
     Computes Obara-Saika factors for electron repulsion integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compFactorsForElectronRepulsion(      CMemBlock2D<double>& osFactors,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair);
    
    /**
     Computes coordinates of combined Gaussian functions, which is obtained by
     applying Gaussian product rule to thee Gaussian functions.
     
     @param wCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimension of Obara-Saika factors vectors.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compCoordinatesForW(      CMemBlock2D<double>& wCoordinates,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nFactors,
                             const CGtoPairsBlock&      braGtoPairsBlock,
                             const CGtoPairsBlock&      ketGtoPairsBlock,
                             const int32_t              nKetPrimPairs,
                             const int32_t              iContrPair);
    
   /**
    Computes vector of distances between center W of combined primitive GTOs
    and combined center P of primitive GTOs pair on bra side.
    
    @param wpDistances the vector of Cartesian R(WP) = W - P distances.
    @param wCoordinates the vector of coordinates for combined Gaussian
           functions.
    @param braGtoPairsBlock the GTOs pairs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compDistancesWP(      CMemBlock2D<double>& wpDistances,
                         const CMemBlock2D<double>& wCoordinates,
                         const CGtoPairsBlock&      braGtoPairsBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              nKetPrimPairs,
                         const int32_t              iContrPair);
    
    /**
     Computes vector of distances between center W of combined primitive GTOs
     and combined center Q of primitive GTOs pair on ket side.
     
     @param wqDistances the vector of Cartesian R(WQ) = W - Q distances.
     @param wCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compDistancesWQ(      CMemBlock2D<double>& wqDistances,
                         const CMemBlock2D<double>& wCoordinates,
                         const CGtoPairsBlock&      braGtoPairsBlock,
                         const CGtoPairsBlock&      ketGtoPairsBlock,
                         const int32_t              nKetPrimPairs,
                         const int32_t              iContrPair);
    
    /**
     Computes vector of effective distances between contracted GTOs pairs center
     center on bra side and combined GTOs pairs centers on ket side.

     @param pqDistances the vector of effective distances.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void compEffectiveDistancesPQ(      CMemBlock<double>& pqDistances,
                                  const CGtoPairsBlock&    braGtoPairsBlock,
                                  const CGtoPairsBlock&    ketGtoPairsBlock,
                                  const bool               isBraEqualKet,
                                  const int32_t            iContrPair);
    
} // intsfunc namespace

#endif /* TwoIntsFunc_hpp */
