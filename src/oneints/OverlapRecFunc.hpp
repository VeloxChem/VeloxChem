//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include<cstdint>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

    /**
    Computes batch of primitive (S||S) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForSS(      CMemBlock2D<double>& primBuffer,
                                CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& abDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (S||P) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForSP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (S||D) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForSD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (S||F) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForSF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (S||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForSG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (P||S) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForPS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& paDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (P||P) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForPP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (P||D) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForPD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (P||F) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForPF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (P||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForPG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (D||S) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForDS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (D||P) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForDP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (D||D) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForDD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (D||F) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForDF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (D||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForDG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (F||S) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (F||P) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (F||D) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (F||F) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (G||S) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForGS(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (G||P) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForGP(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (G||D) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForGD(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (G||F) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForGF(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes batch of primitive (G||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForGG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);


} // ovlrecfunc namespace

