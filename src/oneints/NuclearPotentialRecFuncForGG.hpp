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

namespace npotrecfunc { // npotrecfunc namespace

    /**
    Computes batch of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

    /**
    Computes sub-batch (0,1) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_0_1(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (1,2) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_1_2(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (2,3) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_2_3(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (3,4) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_3_4(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (4,5) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_4_5(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (5,6) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_5_6(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (6,7) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_6_7(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (7,8) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_7_8(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (8,9) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_8_9(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (9,10) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_9_10(      CMemBlock2D<double>& primBuffer,
                                        const CMemBlock2D<double>& auxBuffer,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& paDistances,
                                        const CMemBlock2D<double>& pbDistances,
                                        const CMemBlock2D<double>& pcDistances,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);

    /**
    Computes sub-batch (10,11) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_10_11(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (11,12) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_11_12(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (12,13) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_12_13(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (13,14) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_13_14(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (14,15) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_14_15(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (15,16) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_15_16(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (16,17) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_16_17(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (17,18) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_17_18(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (18,19) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_18_19(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (19,20) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_19_20(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (20,21) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_20_21(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (21,22) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_21_22(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (22,23) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_22_23(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (23,24) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_23_24(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (24,25) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_24_25(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (25,26) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_25_26(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (26,27) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_26_27(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (27,28) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_27_28(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (28,29) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_28_29(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (29,30) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_29_30(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (30,31) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_30_31(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (31,32) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_31_32(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (32,33) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_32_33(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (33,34) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_33_34(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (34,35) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_34_35(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (35,36) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_35_36(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (36,37) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_36_37(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (37,38) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_37_38(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (38,39) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_38_39(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (39,40) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_39_40(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (40,41) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_40_41(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (41,42) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_41_42(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (42,43) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_42_43(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (43,44) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_43_44(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (44,45) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_44_45(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (45,46) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_45_46(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (46,47) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_46_47(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (47,48) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_47_48(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (48,49) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_48_49(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (49,50) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_49_50(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (50,51) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_50_51(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (51,52) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_51_52(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (52,53) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_52_53(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (53,54) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_53_54(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (54,55) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_54_55(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (55,56) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_55_56(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (56,57) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_56_57(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (57,58) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_57_58(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (58,59) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_58_59(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (59,60) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_59_60(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (60,61) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_60_61(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (61,62) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_61_62(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (62,63) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_62_63(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (63,64) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_63_64(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (64,65) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_64_65(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (65,66) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_65_66(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (66,67) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_66_67(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (67,68) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_67_68(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (68,69) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_68_69(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (69,70) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_69_70(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (70,71) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_70_71(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (71,72) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_71_72(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (72,73) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_72_73(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (73,74) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_73_74(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (74,75) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_74_75(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (75,76) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_75_76(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (76,77) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_76_77(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (77,78) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_77_78(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (78,79) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_78_79(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (79,80) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_79_80(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (80,81) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_80_81(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (81,82) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_81_82(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (82,83) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_82_83(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (83,84) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_83_84(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (84,85) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_84_85(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (85,86) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_85_86(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (86,87) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_86_87(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (87,88) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_87_88(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (88,89) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_88_89(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (89,90) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_89_90(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (90,91) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_90_91(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (91,92) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_91_92(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (92,93) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_92_93(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (93,94) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_93_94(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (94,95) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_94_95(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (95,96) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_95_96(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (96,97) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_96_97(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (97,98) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_97_98(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (98,99) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_98_99(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (99,100) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_99_100(      CMemBlock2D<double>& primBuffer,
                                          const CMemBlock2D<double>& auxBuffer,
                                          const CMemBlock2D<double>& osFactors,
                                          const CMemBlock2D<double>& paDistances,
                                          const CMemBlock2D<double>& pbDistances,
                                          const CMemBlock2D<double>& pcDistances,
                                          const CGtoBlock&           braGtoBlock,
                                          const CGtoBlock&           ketGtoBlock,
                                          const int32_t              iContrGto);

    /**
    Computes sub-batch (100,101) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_100_101(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (101,102) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_101_102(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (102,103) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_102_103(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (103,104) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_103_104(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (104,105) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_104_105(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (105,106) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_105_106(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (106,107) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_106_107(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (107,108) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_107_108(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (108,109) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_108_109(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (109,110) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_109_110(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (110,111) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_110_111(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (111,112) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_111_112(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (112,113) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_112_113(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (113,114) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_113_114(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (114,115) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_114_115(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (115,116) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_115_116(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (116,117) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_116_117(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (117,118) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_117_118(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (118,119) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_118_119(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (119,120) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_119_120(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (120,121) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_120_121(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (121,122) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_121_122(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (122,123) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_122_123(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (123,124) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_123_124(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (124,125) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_124_125(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (125,126) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_125_126(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (126,127) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_126_127(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (127,128) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_127_128(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (128,129) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_128_129(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (129,130) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_129_130(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (130,131) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_130_131(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (131,132) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_131_132(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (132,133) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_132_133(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (133,134) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_133_134(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (134,135) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_134_135(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (135,136) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_135_136(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (136,137) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_136_137(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (137,138) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_137_138(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (138,139) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_138_139(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (139,140) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_139_140(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (140,141) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_140_141(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (141,142) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_141_142(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (142,143) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_142_143(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (143,144) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_143_144(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (144,145) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_144_145(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (145,146) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_145_146(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (146,147) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_146_147(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (147,148) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_147_148(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (148,149) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_148_149(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (149,150) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_149_150(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (150,151) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_150_151(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (151,152) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_151_152(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (152,153) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_152_153(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (153,154) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_153_154(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (154,155) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_154_155(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (155,156) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_155_156(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (156,157) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_156_157(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (157,158) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_157_158(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (158,159) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_158_159(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (159,160) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_159_160(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (160,161) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_160_161(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (161,162) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_161_162(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (162,163) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_162_163(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (163,164) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_163_164(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (164,165) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_164_165(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (165,166) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_165_166(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (166,167) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_166_167(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (167,168) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_167_168(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (168,169) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_168_169(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (169,170) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_169_170(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (170,171) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_170_171(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (171,172) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_171_172(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (172,173) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_172_173(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (173,174) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_173_174(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (174,175) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_174_175(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (175,176) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_175_176(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (176,177) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_176_177(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (177,178) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_177_178(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (178,179) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_178_179(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (179,180) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_179_180(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (180,181) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_180_181(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (181,182) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_181_182(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (182,183) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_182_183(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (183,184) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_183_184(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (184,185) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_184_185(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (185,186) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_185_186(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (186,187) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_186_187(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (187,188) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_187_188(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (188,189) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_188_189(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (189,190) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_189_190(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (190,191) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_190_191(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (191,192) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_191_192(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (192,193) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_192_193(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (193,194) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_193_194(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (194,195) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_194_195(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (195,196) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_195_196(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (196,197) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_196_197(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (197,198) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_197_198(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (198,199) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_198_199(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (199,200) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_199_200(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (200,201) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_200_201(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (201,202) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_201_202(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (202,203) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_202_203(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (203,204) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_203_204(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (204,205) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_204_205(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (205,206) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_205_206(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (206,207) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_206_207(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (207,208) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_207_208(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (208,209) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_208_209(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (209,210) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_209_210(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (210,211) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_210_211(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (211,212) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_211_212(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (212,213) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_212_213(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (213,214) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_213_214(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (214,215) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_214_215(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (215,216) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_215_216(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (216,217) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_216_217(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (217,218) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_217_218(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (218,219) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_218_219(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (219,220) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_219_220(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (220,221) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_220_221(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (221,222) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_221_222(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (222,223) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_222_223(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (223,224) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_223_224(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (224,225) of primitive (G|A|G) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGG_224_225(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);


} // npotrecfunc namespace

