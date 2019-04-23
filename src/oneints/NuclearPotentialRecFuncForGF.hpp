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
    Computes batch of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

    /**
    Computes sub-batch (0,3) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_0_3(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (3,6) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_3_6(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (6,9) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_6_9(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (9,12) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_9_12(      CMemBlock2D<double>& primBuffer,
                                        const CMemBlock2D<double>& auxBuffer,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& paDistances,
                                        const CMemBlock2D<double>& pbDistances,
                                        const CMemBlock2D<double>& pcDistances,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);

    /**
    Computes sub-batch (12,15) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_12_15(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (15,18) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_15_18(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (18,21) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_18_21(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (21,24) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_21_24(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (24,27) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_24_27(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (27,30) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_27_30(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (30,33) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_30_33(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (33,36) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_33_36(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (36,39) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_36_39(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (39,42) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_39_42(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (42,45) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_42_45(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (45,48) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_45_48(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (48,51) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_48_51(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (51,54) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_51_54(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (54,57) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_54_57(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (57,60) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_57_60(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (60,63) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_60_63(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (63,66) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_63_66(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (66,69) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_66_69(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (69,72) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_69_72(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (72,75) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_72_75(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (75,78) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_75_78(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (78,81) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_78_81(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (81,84) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_81_84(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (84,87) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_84_87(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (87,90) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_87_90(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (90,93) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_90_93(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (93,96) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_93_96(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (96,99) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_96_99(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (99,102) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_99_102(      CMemBlock2D<double>& primBuffer,
                                          const CMemBlock2D<double>& auxBuffer,
                                          const CMemBlock2D<double>& osFactors,
                                          const CMemBlock2D<double>& paDistances,
                                          const CMemBlock2D<double>& pbDistances,
                                          const CMemBlock2D<double>& pcDistances,
                                          const CGtoBlock&           braGtoBlock,
                                          const CGtoBlock&           ketGtoBlock,
                                          const int32_t              iContrGto);

    /**
    Computes sub-batch (102,105) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_102_105(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (105,108) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_105_108(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (108,111) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_108_111(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (111,114) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_111_114(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (114,117) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_114_117(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (117,120) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_117_120(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (120,123) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_120_123(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (123,126) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_123_126(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (126,129) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_126_129(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (129,132) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_129_132(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (132,135) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_132_135(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (135,138) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_135_138(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (138,141) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_138_141(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (141,144) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_141_144(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (144,147) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_144_147(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);

    /**
    Computes sub-batch (147,150) of primitive (G|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForGF_147_150(      CMemBlock2D<double>& primBuffer,
                                           const CMemBlock2D<double>& auxBuffer,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& paDistances,
                                           const CMemBlock2D<double>& pbDistances,
                                           const CMemBlock2D<double>& pcDistances,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);


} // npotrecfunc namespace

