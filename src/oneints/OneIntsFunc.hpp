//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OneIntsFunc_hpp
#define OneIntsFunc_hpp

#include <cstdint>
#include <cmath>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

namespace intsfunc { // intsfunc namespace
    
    /**
     Computes distances between specific contracted GTO on bra side with all
     contracted GTOs on ket side.

     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesAB(      CMemBlock2D<double>& abDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for overlap integrals.

     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForOverlap(      CMemBlock2D<double>& osFactors,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center A of primitive GTO on bra side.
     
     @param paDistances the vector of Cartesian R(PA) = P - A distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPA(      CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center B of primitive GTO on ket side.
     
     @param pbDistances the vector of Cartesian R(PB) = P - B distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPB(      CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
} // intsfunc namespace

#endif /* OneIntsFunc_hpp */
