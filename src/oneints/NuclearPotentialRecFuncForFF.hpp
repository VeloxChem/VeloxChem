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
    Computes batch of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

    /**
    Computes sub-batch (0,3) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_0_3(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (3,6) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_3_6(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (6,9) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_6_9(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (9,12) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_9_12(      CMemBlock2D<double>& primBuffer,
                                        const CMemBlock2D<double>& auxBuffer,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& paDistances,
                                        const CMemBlock2D<double>& pbDistances,
                                        const CMemBlock2D<double>& pcDistances,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);

    /**
    Computes sub-batch (12,15) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_12_15(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (15,18) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_15_18(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (18,21) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_18_21(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (21,24) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_21_24(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (24,27) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_24_27(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (27,30) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_27_30(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (30,33) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_30_33(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (33,36) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_33_36(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (36,39) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_36_39(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (39,42) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_39_42(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (42,45) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_42_45(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (45,48) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_45_48(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (48,51) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_48_51(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (51,54) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_51_54(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (54,57) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_54_57(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (57,60) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_57_60(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (60,63) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_60_63(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (63,66) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_63_66(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (66,69) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_66_69(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (69,72) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_69_72(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (72,75) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_72_75(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (75,78) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_75_78(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (78,81) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_78_81(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (81,84) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_81_84(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (84,87) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_84_87(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (87,90) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_87_90(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (90,93) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_90_93(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (93,96) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_93_96(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (96,98) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_96_98(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (98,100) of primitive (F|A|F) nuclear potential integrals and stores
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
    void compNuclearPotentialForFF_98_100(      CMemBlock2D<double>& primBuffer,
                                          const CMemBlock2D<double>& auxBuffer,
                                          const CMemBlock2D<double>& osFactors,
                                          const CMemBlock2D<double>& paDistances,
                                          const CMemBlock2D<double>& pbDistances,
                                          const CMemBlock2D<double>& pcDistances,
                                          const CGtoBlock&           braGtoBlock,
                                          const CGtoBlock&           ketGtoBlock,
                                          const int32_t              iContrGto);


} // npotrecfunc namespace

