//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OverlapRecFunc_hpp
#define OverlapRecFunc_hpp

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace
    
    /**
     Computes batch of primitive (S|S) overlap integrals and stores results in
     primitives buffer.

     @param primBuffer the primitives buffer.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForSS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& abDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);
    
    /**
     Computes batch of primitive (S|P) overlap integrals and stores results in
     primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForSP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);
    
    /**
     Computes batch of primitive (S|D) overlap integrals and stores results in
     primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForSD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);
    
    /**
     Computes batch of primitive (S|F) overlap integrals and stores results in
     primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForSF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);
    
    /**
     Computes batch of primitive (S|G) overlap integrals and stores results in
     primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForSG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);
    
    /**
     Computes batch of primitive (P|S) overlap integrals and stores results in
     primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compOverlapForPS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& paDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

} // ovlrecfunc namespace

#endif /* OverlapRecFunc_hpp */
